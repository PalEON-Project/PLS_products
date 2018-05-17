library(raster)
library(dplyr)
library(ncdf4)

if(shared_params_in_cell)
    load('point_with_biomass_shared.Rda') else load('point_with_biomass.Rda')

if(!'cell' %in% names(mw)) 
    mw <- mw %>% add_cells_to_dataset() 

xy <- xyFromCell(base_raster, getValues(base_raster))
grid <- tibble(x = xy[,1], y = xy[,2], cell = getValues(base_raster))

mask <- nc_open(file.path(conversions_data_dir, 'paleonmask.nc'))
regions <- ncvar_get(mask, 'subregion', c(1,1),c(-1,-1))

x <- matrix(ncvar_get(mask, 'x', c(1),c(-1)), nrow = nrow(regions),
            ncol = ncol(regions))
y <- matrix(ncvar_get(mask, 'y', c(1),c(-1)),nrow = nrow(regions),
            ncol = ncol(regions), byrow = TRUE) 

west <- c(regions %in% paleon_regions_west)
paleon <- c(regions %in% paleon_regions)
pred_grid <- data.frame(x = c(x), y = c(y))
pred_grid_paleon <- pred_grid[paleon, ]
pred_grid_west <- pred_grid[west, ]


## total biomass

k_pot <- 2500
k_occ <- 1500  ## should probably be larger but takes a long time to fit; explore this further

## compute average biomass per tree; this uses two-tree points with one tree with missing diameter
biomass_avg <- mw %>% dplyr::select(biomass1, biomass2) %>% as.matrix(.)  %>% apply(1, mean, na.rm = TRUE)

mw <- mw %>% mutate(biomass = ifelse(num_trees == 0, 0,
                                        ifelse(num_trees == 1, biomass1, biomass_avg)))

cell_full <- mw %>% filter(!(is.na(biomass) | is.na(density_for_biomass))) %>% group_by(cell) %>% summarize(total = n())

cell_occ <- mw %>% filter(!(is.na(biomass) | is.na(density_for_biomass)) & biomass > 0) %>% 
    group_by(cell) %>%
    summarize(avg = mean(biomass*density_for_biomass/kg_per_Mg),
              geom_avg = mean(log(biomass*density_for_biomass/kg_per_Mg)),
              count = n())

## note: given we fit stat model for potential biomass on log scale, can't really scale
## by count unless we use geometric avg

cell_full <- cell_full %>% left_join(cell_occ, by = c("cell" = "cell")) %>% left_join(grid, by = c("cell" = "cell")) %>% mutate(count = ifelse(is.na(count), 0 , count))

## explore bam vs gam further
test <- fit(cell_full, newdata = pred_grid_west, k_occ = k_occ, k_pot = k_pot, use_bam = TRUE)

## comparisons

if(F) {
par(mfrow=c(2,2))
raw <- mw %>% filter(!(is.na(biomass) | is.na(density_for_biomass))) %>% group_by(cell) %>%
    summarize(biom = mean(biomass*density_for_biomass / kg_per_Mg, na.rm = TRUE))
raw <- raw %>% right_join(grid, by = c("cell" = "cell")) %>% mutate(biom_tr = ifelse(biom > 300, 300, biom))
image.plot(sort(unique(raw$x)), sort(unique(raw$y)), matrix(raw$biom_tr, 296, 180)[,180:1],xlim=c(-70000,1100000))

pmap2(cell_full$avg, adf(cell_full[,c('x','y')]), cex=.5,zlim=c(0,300))
pmap2(exp(cell_full$geom_avg), adf(cell_full[,c('x','y')]), cex=0.5,zlim=c(0,300))



scaling=8000
gamma=1
data=cell_full
newdata = pred_grid_west
newdata[ , c('x','y')] <- newdata[ , c('x','y')] / scaling
    data[ , c('x','y')] <- data[ , c('x','y')] / scaling 
     data$z <- cbind(data$count, data$total - data$count)

        ## default method is "fREML" which doesn't seem to respond to 'gamma'
        ## and is not the same method as used in gam()
       system.time( occ <- bam(z ~ s(x, y, k = 1500), data = data, family = 'binomial',
                   gamma = gamma, method = "GCV.Cp"))       
pr_occ <- predict(occ, newdata = newdata, type = 'response')
  data <- data %>% filter(count > 0)
   
        data$z <- data$avg
    pot_aa <- gam(z ~ s(x,y, k = 2500), data = data, weights = count,
               gamma = gamma, method = "GCV.Cp")
        data$z <- log(data$avg)  
    pot_a <- gam(z ~ s(x,y, k = 2500), data = data, weights = count,
               gamma = gamma, method = "GCV.Cp")
      data$z <- data$geom_avg
    pot_g <- gam(z ~ s(x,y, k = 2500), data = data, weights = count,
               gamma = gamma, method = "GCV.Cp")
    pr_pot_aa <- (predict(pot_aa, newdata= newdata, type='response'))
    pr_pot_a <- exp(predict(pot_a, newdata= newdata, type='response'))
    pr_pot_g <- exp(predict(pot_g, newdata= newdata, type='response'))

b_aa = pr_occ * pr_pot_aa
b_a = pr_occ * pr_pot_a
b_g = pr_occ * pr_pot_g

par(mfrow=c(1,3))
pmap2(b_aa, newdata[,c('x','y')], cex = .5,zlim=c(0,300))
pmap2(b_a, newdata[,c('x','y')], cex = .5,zlim=c(0,300))
pmap2(b_g, newdata[,c('x','y')], cex = .5,zlim=c(0,300))

}

## CV

set.seed(1)
cells <- sample(unique(cell_full$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))

gamma_values <- c(.01, .03, .1, .3, .7, 1, 2, 3, 10, 30, 100)

cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

pred_arith_fit <- pred_geom_fit <- se_arith_fit <- se_geom_fit <- matrix(0, nrow(cell_full), length(gamma_values))
for(i in 1:n_folds) {
    train <- cell_full %>% filter(fold != i)
    test <- cell_full %>% filter(fold == i)

    for(j in seq_along(gamma_values)) {
        tmp <- fit(train, newdata = test, gamma = gamma_values[j], k_occ = k_occ, k_pot = k_pot, use_bam = TRUE)
        pred_arith_fit[cell_full$fold == i, j] <- tmp$mean
        se_arith_fit[cell_full$fold == i, j] <- tmp$sd
        
        tmp <- fit(train, newdata = test, do_log = FALSE, gamma = gamma_values[j], k_occ = k_occ, k_pot = k_pot, use_bam = TRUE)
        pred_geom_fit[cell_full$fold == i, j] <- tmp$mean
        se_geom_fit[cell_full$fold == i, j] <- tmp$sd
    }
}

wgt_sse <- function(n, y, yhat) {
    sum(n * (y - yhat)^2)
}

# geom fit looks a bit better and gamma=1 or slightly higher seems best
apply(pred_arith_fit, 2, function(x) wgt_sse(cell_full$total, cell_full$avg, x))
apply(pred_geom_fit, 2, function(x) wgt_sse(cell_full$total, cell_full$avg, x))

## biomass by taxon

k_occ <- 1500
k_pot <- 2500

taxon = "Oak"

calc_biomass_taxon <- function(num_trees, biomass1, biomass2, density, L3_tree1, L3_tree2, taxon) {
    biomass <- rep(0, length(num_trees))
    
    cond <- num_trees == 1 & L3_tree1 == taxon
    biomass[cond] <- biomass1[cond] * density[cond] / kg_per_Mg
    ## assign half biomass to taxon if either tree is the taxon
    cond <- num_trees == 2 & L3_tree1 == taxon 
    biomass[cond] <- biomass1[cond] * density[cond]  / (kg_per_Mg * 2)  ## 2 to account for taxon represents half the density
    ## if two trees of same taxon, the addition should handle this
    cond <- num_trees == 2 & L3_tree2 == taxon 
    biomass[cond] <- biomass[cond] + biomass2[cond] * density[cond]  / (kg_per_Mg * 2)

    ## handle case of two trees same taxon but one biomass is missing; use single biomass as the per-tree estimate
    cond <- (num_trees == 2 & L3_tree1 == taxon & L3_tree2 == taxon & is.na(biomass1)) |
        (num_trees == 2 & L3_tree1 == taxon & L3_tree2 == taxon & is.na(biomass2))
    biomass[cond] <- biomass2[cond] * density[cond]  / kg_per_Mg
    
    return(biomass)
}

tmp <- mw %>% mutate(biomass_focal = calc_biomass_taxon(num_trees, biomass1, biomass2, density_for_biomass, L3_tree1, L3_tree2, taxon))

cell_full <- tmp %>% filter(!(is.na(biomass_focal))) %>% group_by(cell) %>% summarize(total = n())

cell_occ <- tmp %>% filter(!is.na(biomass_focal) & biomass_focal > 0) %>% 
    group_by(cell) %>%
    summarize(avg = mean(biomass_focal),
              geom_avg = mean(log(biomass_focal)),
              count = n())

cell_full <- cell_full %>% left_join(cell_occ, by = c("cell" = "cell")) %>% left_join(grid, by = c("cell" = "cell")) %>% mutate(count = ifelse(is.na(count), 0 , count))

test <- fit(cell_full, newdata = pred_grid_west, k_occ = k_occ, k_pot = k_pot, use_bam = TRUE)

## CV

set.seed(1)
cells <- sample(unique(cell_full$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))

gamma_values <- c(.01, .03, .1, .3, .7, 1, 2, 3, 10, 30, 100)

cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

pred_arith_fit <- pred_geom_fit <- se_arith_fit <- se_geom_fit <- matrix(0, nrow(cell_full), length(gamma_values))
for(i in 1:n_folds) {
    train <- cell_full %>% filter(fold != i)
    test <- cell_full %>% filter(fold == i)

    for(j in seq_along(gamma_values)) {
        tmp <- fit(train, newdata = test, gamma = gamma_values[j], k_occ = k_occ, k_pot = k_pot, use_bam = TRUE)
        pred_arith_fit[cell_full$fold == i, j] <- tmp$mean
        se_arith_fit[cell_full$fold == i, j] <- tmp$sd
        
        tmp <- fit(train, newdata = test, do_log = FALSE, gamma = gamma_values[j], k_occ = k_occ, k_pot = k_pot, use_bam = TRUE)
        pred_geom_fit[cell_full$fold == i, j] <- tmp$mean
        se_geom_fit[cell_full$fold == i, j] <- tmp$sd
    }
}

wgt_sse <- function(n, y, yhat) {
    sum(n * (y - yhat)^2)
}

# geom fit looks a bit better and gamma=1 or slightly higher seems best
apply(pred_arith_fit, 2, function(x) wgt_sse(cell_full$total, cell_full$avg, x))
apply(pred_geom_fit, 2, function(x) wgt_sse(cell_full$total, cell_full$avg, x))
