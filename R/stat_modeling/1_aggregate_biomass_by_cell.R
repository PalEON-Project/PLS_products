library(raster)
library(dplyr)
library(ncdf4)

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

test <- fit(cell_full, newdata = pred_grid_west, k_occ = k_occ, k_pot = k_pot, use_bam = TRUE)

## c-v

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

k_occ <- 3000
k_pot <- 1000

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
    cond <- num_trees == 2 & L3_tree1 == taxon & L3_tree2 == taxon & is.na(biomass1)
    biomass[cond] <- biomass2[cond] * density[cond]  / kg_per_Mg
    cond <- num_trees == 2 & L3_tree1 == taxon & L3_tree2 == taxon & is.na(biomass2)
    biomass[cond] <- biomass1[cond] * density[cond]  / kg_per_Mg
    
    return(dens)
}

## HERE

mw <- mw %>% summarize(biomass = ifelse(num_trees == 0, 0,
                                        ifelse(num_trees == 1, biomass1,
                                               apply(as.matrix(mw$biomass1, mw$biomass2), 1,
                                                     mean, na.rm = TRUE)))

cell_full <- mw %>% filter(!(is.na(biomass) | is.na(density_for_biomass))) %>% group_by(cell) %>% summarize(total = n())

cell_occ <- mw %>% filter(!is.na(biomass) & !is.na(density_for_biomass) & biomass > 0) %>% 
    group_by(cell) %>%
    summarize(avg = mean(biomass*density_for_biomass/kg_per_Mg),
              geom_avg = mean(log(biomass*density_for_biomass/kg_per_Mg)),
              count = n())


## OLD
# KH adding the code that has been used in the past to aggregate biomass, density, etc up to the grid cell:

spec.table <- read.csv(paste0('outputs/biomass_no_na_pointwise.ests_inilmi','_v',version, '.csv'))

#-------------------------Paleon gridding----------------------------------------------

# These are not the full tables since the include only the Paleon grid cells with points in the database.
#dcast rearranges the spec.table data by x, y and cell
count.table <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm=TRUE, value.var = 'count')

count.table.pt <- dcast(spec.table, Pointx + Pointy + cell ~ spec, sum, na.rm=TRUE, value.var = 'count')

unique.len <- function(x){length(unique(x))}

#melt data by by x, y, cell as rows and spec as columns and provide the sum or the number of unique points
biomass.trees  <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm=TRUE, value.var = 'count')

#biomass.trees represnts the number of trees in each category

biomass.points <- dcast(spec.table, x + y + cell ~ spec, unique.len, value.var = 'point')
#biomass.points represents the number of unique points
# This is to get the points with trees:
spec.adj <- spec.table
spec.adj$tree <- spec.table$spec %in% "No tree"

treed.points <- dcast(spec.adj, x + y + cell ~ tree, unique.len, value.var = 'point')
treed.points <- treed.points[order(treed.points$cell),]
colnames(treed.points) <- c('x', 'y', 'cell', 'tree', 'no tree')

# Now get the total number of plots per cell:
plots_per_cell <- dcast(spec.table, x + y + cell ~ ., unique.len, value.var = 'point')
plots_per_cell <- plots_per_cell[order(plots_per_cell$cell),]

points_per_cell_df <- data.frame(plots_per_cell,
                                 treed = treed.points$tree,
                                 untreed = treed.points$`no tree`)

colnames(points_per_cell_df)[4:5] <- c("total_pls_points", "total_treed_points")


# points by cell is the # of points 
points.by.cell <- rowSums(biomass.points[, 4:ncol(biomass.points)], na.rm=TRUE)
trees.by.cell  <- rowSums(biomass.trees[,!colnames(biomass.trees) %in% c('x', 'y', 'cell', 'No tree','Wet', 'Water')], na.rm=TRUE)

#calculate the sum of total density, basal area, biomass & diameter by cell and species
density.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'density')
basal.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'basal')
biomass.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'biom')
diam.table <-  dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'diams')



#calculate standard deviations of density, basal area, biomass, and diameters
density.sd.table <- dcast(spec.table, x + y  + cell ~ spec, sd, na.rm=TRUE, value.var = 'density')
basal.sd.table <- dcast(spec.table, x + y  + cell ~ spec, sd, na.rm=TRUE, value.var = 'basal')
diam.sd.table <-  dcast(spec.table, x + y  + cell ~ spec, sd, na.rm=TRUE, value.var = 'diams')
biom.sd.table <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm=TRUE, value.var = 'biom')

#  The function averages the estimates a single value for density, basal area, and biomass based on the points.by.cell
normalize <- function(x, mult = 2, value = points.by.cell) {
  x[,4:ncol(x)] <-  x[,4:ncol(x)] / value *mult; x}

density.table <- normalize(density.table)
basal.table <- normalize(basal.table)
diam.table <- normalize(diam.table, mult = 2.54, trees.by.cell)
biomass.table <- normalize(biomass.table)

# output preliminary density table to look at:
write.csv(density.table, "outputs/density.table_test.csv")

# get a density total value
density.table$total = rowSums(density.table[,4:ncol(density.table)], na.rm=TRUE)

mass.table <- merge(biomass.table, points_per_cell_df[,c("x", "y", "cell", "total_pls_points", "total_treed_points")], by = c("x", "y", "cell"))

#  We want rasterized versions of these tables with sums:
rast.fun <- function(x){
  
 to_grid <- data.frame(cell = x$cell, 
                  total = rowSums(x[,4:ncol(x)], na.rm=TRUE))
  
  empty <- rep(NA, ncell(base.rast))
  empty[to_grid$cell] <- to_grid$total
  setValues(base.rast, empty)
  }

base.rast <- raster(xmn = -71000, xmx = 2297000, ncols=296,
                    ymn = 58000,  ymx = 1498000, nrows = 180,
                    crs = '+init=epsg:3175')

#create total rasters for densi, basal area, diameters, biomass
dens     <- rast.fun(density.table)
count.up <- rast.fun(count.table)
basal    <- rast.fun(basal.table)
mdiam    <- rast.fun(diam.table); mdiam[mdiam==0] <- NA
biomass <- rast.fun(biomass.table)

writeRaster(biomass, "data/biomass.grd", overwrite=TRUE)
writeRaster(dens, "data/dens.grd", overwrite=TRUE)


#to get sd for these
rowset <- cbind(1:(nrow(spec.table)/2), (nrow(spec.table)/2+1):nrow(spec.table))

sd.table <- data.frame(cell = spec.table$cell[rowset[,1]],
                       density = spec.table$density[rowset[,1]] * 2,
                       basal   = spec.table$basal[rowset[,1]] + spec.table$basal[rowset[,2]],
                       biomass = spec.table$biom[rowset[,1]] + spec.table$biom[rowset[,2]],
                       diam = spec.table$diam[rowset[,1]] + spec.table$diam[rowset[,2]])

dens.sd <- dcast(sd.table, cell ~ ., fun.aggregate=sd, na.rm=TRUE, value.var = 'density')
basal.sd <- dcast(sd.table, cell ~ ., sd, na.rm=TRUE, value.var = 'basal')
biomass.sd <- dcast(sd.table, cell ~ ., sd, na.rm=TRUE, value.var = 'biomass')
diam.sd <- dcast(sd.table, cell ~., fun.aggregate = sd, na.rm = TRUE, value.var = 'diam')
colnames(dens.sd) <-c("cell", "dens.sd")
colnames(basal.sd) <-c("cell", "basal.sd")
colnames(biomass.sd) <-c("cell", "biomass.sd")
colnames(diam.sd) <-c("cell", "diam.sd")

#combine sd's with the tables
density.table <- merge(density.table, dens.sd, by = "cell")
diam.table <- merge(diam.table, diam.sd, by = "cell")



data.table <- data.frame(xyFromCell(dens, 1:ncell(dens)), 
                         stem.dens = getValues(dens),
                         basal.area = getValues(basal),
                         biomass = getValues(biomass),
                         diam = getValues(mdiam))

data.table <- data.table[!is.na(data.table[,3]), ] # remove places missing cells
write.csv(data.table, paste0('data/density.basal.biomass_alb','_v',1, '.csv'))


#  Now we need to add zero cells to the dataframe:
comp.table <- data.frame(xyFromCell(base.rast, 1:ncell(base.rast)), 
                         cell = 1:ncell(base.rast),
                         matrix(ncol = ncol(biomass.table)-3, nrow = ncell(base.rast)))

colnames(comp.table)[4:ncol(biomass.table)] <- colnames(biomass.table)[4:ncol(comp.table)]

reform <- function(x){
  comp.table[x$cell,] <- x
  comp.table
}


biomass.full <- reform(biomass.table)
density.full <- reform(density.table)
diameter.full <- reform(diam.table)
composition.table <- reform(basal.table)


#make composition table
composition.table <- composition.table[,!names(composition.table) %in% c('Water', "Wet", "No tree", "plss_pts")]
composition.table[,4:ncol(composition.table)] <- composition.table[,4:ncol(composition.table)]/rowSums(composition.table[,4:ncol(composition.table)], na.rm=TRUE)

# composition table without paleon grid NA values:
comp.inil <- composition.table[complete.cases(composition.table),]


#extract biomass.full table by the extent to get a datatable


biomass.df.test <- extract(biomass, extent(ind_il), cellnumbers=TRUE, xy=TRUE)
#biomass.df.test has the cell numbers of the extent of indiana, so I will use these

#get the biomass.full cells out of biomass.full
biomass.df.test.cells <- data.frame(biomass.df.test[,1])
colnames(biomass.df.test.cells) <- 'cell'
biomass.indiana <- merge(biomass.full, biomass.df.test.cells, by ="cell")

#to get biomass.pts and tree.pts, use merge(x, y, by, all=TRUE)
cells.xy.ind <- biomass.indiana[,1:3] #xy and cells in indiana
biomass.points.ind <- merge(cells.xy.ind, biomass.points, all=TRUE)
count.trees.ind <- merge(cells.xy.ind, count.table, all=TRUE)


add.v <- function(x, name, row.names){
  
  #  Quick file name formatter:
  
    p.ext <- '_alb'
  
  write.csv(x, paste0('data/outputs/', name,  '_v',version, '.csv'), row.names = row.names)
}

add.v(count.table, 'plss_counts', row.names=FALSE)
add.v(biomass.points, 'plss_points', row.names=FALSE)
add.v(count.table, "plss_composition", row.names=FALSE)
add.v(count.table, "plss_inil_composition", row.names=FALSE)
biomass.table$plsspts_cell <- points.by.cell

add.v(density.table, 'plss_density', row.names=FALSE)
add.v(basal.table, 'plss_basal', row.names=FALSE)
add.v(biomass.table, 'plss_biomass', row.names=FALSE)

add.v(biomass.full,  'plss_spec_biomass', row.names=FALSE) #full species biomass for indiana
add.v(density.full, 'plss_spec_density', row.names=FALSE)
add.v(diameter.full, 'plss_spec_diam', row.names=FALSE)
add.v(diam.table, 'plss_diam', row.names=FALSE)
