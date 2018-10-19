## Fit statistical model to smooth the raw cell-level biomass via cross-validation
## to determine best upper-bound on amount of spatial smoothing.

## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

## note that in cases with one tree with 0 or NA diameter and therefore missing biomass
## we use the biomass for the other tree as the representative value
## this assumes that missingness of diameter is non-informative
## 1384 2-tree points with tree2 diam missing
## 1317 2-tree points with tree1 diam missing

library(dplyr)

load(file.path(interim_results_dir, 'cell_with_biomass_grid.Rda'))

library(doParallel)
if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}
registerDoParallel(cores = n_cores)

biomass_avg <- mw %>% dplyr::select(biomass1, biomass2) %>% as.matrix(.) %>%
    apply(1, mean, na.rm = TRUE)

## add total point-level biomass to dataset
mw <- mw %>% mutate(biomass = ifelse(num_trees == 0, 0,
                              ifelse(num_trees == 1, biomass1, biomass_avg)))

## total points per cell
cell_full <- mw %>% filter(!(is.na(biomass) | is.na(density_for_biomass))) %>%
    group_by(cell) %>% summarize(points_total = n())

## biomass stats averaged over occupied points
cell_occ <- mw %>% filter(!(is.na(biomass) | is.na(density_for_biomass)) & biomass > 0) %>% 
    group_by(cell) %>%
    summarize(avg = mean(biomass*density_for_biomass/kg_per_Mg),
              geom_avg = mean(log(biomass*density_for_biomass/kg_per_Mg)),
              points_occ = n())

## note: given we fit stat model for potential biomass on log scale, can't really scale
## by count_occ unless we use geometric or arithmetic average, not log of arithmetic average
## but in cross-validation log of arithmetic average seems to work better

## should have 'total', 'occupied', and biomass stats (for occupied points)
cell_full <- cell_full %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
    left_join(grid, by = c("cell" = "cell")) %>%
    mutate(points_occ = ifelse(is.na(points_occ), 0, points_occ))

set.seed(1)
cells <- sample(cell_full$cell, replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))

cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

pred_occ <- matrix(0, nrow(cell_full), length(k_occ_cv))
dimnames(pred_occ)[[2]] <- k_occ_cv
pred_pot_arith <- pred_pot_larith <- matrix(0, nrow(cell_full), length(k_pot_cv))
dimnames(pred_pot_arith)[[2]] <- dimnames(pred_pot_larith)[[2]] <- k_pot_cv

n_folds <- max(cell_full$fold)
output <- foreach(i = seq_len(n_folds)) %dopar% {
    train <- cell_full %>% filter(fold != i)
    test <- cell_full %>% filter(fold == i)
    
    po <- fit(train, newdata = test, k_occ = k_occ_cv, unc = FALSE, use_bam = TRUE)
    ppa <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', unc = FALSE, use_bam = TRUE)
    ppl <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', unc = FALSE, use_bam = TRUE)
    cat("n_fold: ", i, " ", date(), "\n")
    list(po, ppa, ppl)
}
for(i in seq_len(n_folds)) {
    pred_occ[cell_full$fold == i, ] <- output[[i]][[1]]$pred_occ
    pred_pot_arith[cell_full$fold == i, ] <- output[[i]][[2]]$pred_pot
    pred_pot_larith[cell_full$fold == i, ] <- output[[i]][[3]]$pred_pot
}

## Assess results

y <- cell_full$avg*cell_full$points_occ/cell_full$points_total ## actual average biomass over all cells
y[is.na(y)] <- 0   # cells with no points with trees (since $avg will be NA)

critArith <- calc_cv_criterion(pred_occ, pred_pot_arith, cell_full$points_total,
                               y, cv_max_biomass)
critLogArith <- calc_cv_criterion(pred_occ, pred_pot_larith, cell_full$points_total,
                                  y, cv_max_biomass)

save(critArith, critLogArith, pred_occ, pred_pot_arith, pred_pot_larith, file = file.path(interim_results_dir,
                                                        'cv_total_biomass.Rda'))
