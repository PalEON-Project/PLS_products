## Fit statistical model to smooth the raw point level density.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.

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

set.seed(1)
cells <- sample(unique(cell_full$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))
    
cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}

library(doParallel)
registerDoParallel(cores = n_cores)

results <- fit_cv_total(cell_full, k_occ_cv, k_pot_cv)

## assess results
y <- cell_full$avg*cell_full$points_occ/cell_full$points_total ## actual average density over all cells

critArith <- calc_cv_criterion(results$pred_occ, results$pred_pot_arith, cell_full$points_total, y, cv_max_density)
critLogArith <- calc_cv_criterion(results$pred_occ, results$pred_pot_larith, cell_full$points_total, y, cv_max_density)

save(critArith, critLogArith, results, file = file.path(interim_results_dir, 'cv_total_density.Rda'))

