## Fit statistical model to smooth the raw point level density.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.

## Run-time (without cross-validation) approximately with k values of occ: 500, pot: 3500

if(!exists('k_pot_total'))
    stop("Must specify 'k_pot_total'")
if(!exists('k_occ_total'))
    stop("Must specify 'k_occ_total'")

## note that in cases with one tree with 0 or NA diameter and therefore missing density
## we use the density for the other tree as the representative value
## this assumes that missingness of diameter is non-informative
## 1384 2-tree points with tree2 diam missing
## 1317 2-tree points with tree1 diam missing

## total points per cell
cell_full <- mw %>% filter(!is.na(density)) %>%
    group_by(cell) %>% summarize(points_total = n())

## density stats averaged over occupied points
cell_occ <- mw %>% filter(!is.na(density) & density > 0) %>% 
    group_by(cell) %>%
    summarize(avg = mean(density),
              geom_avg = mean(log(density)),
              points_occ = n())

## note: given we fit stat model for potential density on log scale, can't really scale
## by count_occ unless we use geometric or arithmetic average, not log of arithmetic average
## but in cross-validation log of arithmetic average seems to work better for biomass
## so assuming that is the case for density for the moment as well

## should have 'total', 'occupied', and density stats (for occupied points)
cell_full <- cell_full %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
    left_join(grid, by = c("cell" = "cell")) %>%
    mutate(points_occ = ifelse(is.na(points_occ), 0, points_occ))

if(do_cv) {

    k_occ <- c(100,250,500,1000,1500,2000,2500)
    k_pot = c(100,250,500,1000,1500,2000,2500,3000,3500)
    
    set.seed(1)
    cells <- sample(unique(cell_full$cell), replace = FALSE)
    folds <- rep(1:n_folds, length.out = length(cells))
    
    cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))
    stop("this CV code won't work yet")
    results <- fit_cv(cell_full, k_occ, k_pot, n_cores)

    ## assess results
    
    y <- cell_full$avg*cell_full$points_occ/cell_full$points_total ## actual average density over all cells
    y[is.na(y)] <- 0
    y[y > mx] <- mx


    critArith <- calc_cv_criterion(results$pred_occ, results$pred_pot_arith, cell_full$points_total,
                                   cell_full$avg*cell_full$count/cell_full$total, 200)
    critLogArith <- calc_cv_criterion(results$pred_occ, results$pred_pot_larith, cell_full$points_total,
                                   cell_full$avg*cell_full$count/cell_full$total, 200)

}

## fit stats model
density_total <- fit(cell_full, newdata = pred_grid_west, k_occ = k_occ_total, k_pot = k_pot_total, return_model = TRUE, unc = TRUE, type_pot = 'log_arith', num_draws = n_stat_samples, save_draws = TRUE)

save(density_total, file = file.path(interim_results_dir, 'fitted_total_density.Rda'))
