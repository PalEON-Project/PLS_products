## Fit statistical model to smooth the raw cell-level density via cross-validation
## to determine best upper-bound on amount of spatial smoothing.


## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.

library(dplyr)

load(file.path(interim_results_dir, 'cell_with_density_grid.Rda'))

library(doParallel)
if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}
registerDoParallel(cores = n_cores)

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
cells <- sample(cell_full$cell, replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))
    
cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

pred_occ <- matrix(0, nrow(cell_full), length(k_occ_cv))
dimnames(pred_occ)[[2]] <- k_occ_cv
pred_pot_arith <- pred_pot_larith <- sig2_arith <- sig2_larith <- matrix(0, nrow(cell_full), length(k_pot_cv))
dimnames(pred_pot_arith)[[2]] <- dimnames(pred_pot_larith)[[2]] <- dimnames(sig2_arith)[[2]] <-
    dimnames(sig2_larith)[[2]] <-k_pot_cv

draws_logocc <- array(0, c(nrow(cell_full), length(k_occ_cv), n_stat_samples))
draws_logpot_arith1 <- draws_logpot_larith1 <- draws_logpot_arith70 <- draws_logpot_larith70 <- array(0, c(nrow(cell_full), length(k_pot_cv), n_stat_samples))


n_folds <- max(cell_full$fold)
output <- foreach(i = seq_len(n_folds)) %dopar% {
    train <- cell_full %>% filter(fold != i)
    test <- cell_full %>% filter(fold == i)
    
    po <- fit(train, newdata = test, k_occ = k_occ_cv, unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
    ppa <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', unc = TRUE, save_draws = TRUE, use_bam = TRUE, num_draws = n_stat_samples)
    ppl <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', unc = TRUE, save_draws = TRUE, use_bam = TRUE, num_draws = n_stat_samples)
    cat("n_fold: ", i, " ", date(), "\n")
    list(po = po, ppa = ppa, ppl = ppl)
}
for(i in seq_len(n_folds)) {
    subn <- sum(cell_full$fold == i)
    pred_occ[cell_full$fold == i, ] <- output[[i]]$po$pred_occ
    pred_pot_arith[cell_full$fold == i, ] <- output[[i]]$ppa$pred_pot
    pred_pot_larith[cell_full$fold == i, ] <- output[[i]]$ppl$pred_pot
    draws_logocc[cell_full$fold == i, , ] <- output[[i]]$po$draws_logocc
    draws_logpot_arith[cell_full$fold == i, , ] <- output[[i]]$ppa$draws_logpot
    draws_logpot_larith[cell_full$fold == i, , ] <- output[[i]]$ppl$draws_logpot
    sig2_arith[cell_full$fold == i, ] <- rep(output[[i]]$ppa$model_pot, each = subn)
    sig2_larith[cell_full$fold == i, ] <- rep(output[[i]]$ppl$model_pot, each = subn)
}


## assess results
y <- cell_full$avg*cell_full$points_occ/cell_full$points_total ## actual average density over all cells
y[is.na(y)] <- 0   # cells with no points with trees (since $avg will be NA)

crit_arith <- calc_cv_criterion(pred_occ, pred_pot_arith, cell_full$points_total,
                               y, cv_max_density)
crit_larith <- calc_cv_criterion(pred_occ, pred_pot_larith, cell_full$points_total,
                                  y, cv_max_density)

cell_full$obs <- y
crit_arith <- c(list(point = crit_arith),
                 calc_cov_criterion(draws_logocc, draws_logpot_arith, sig2 = sig2_arith,
                                    cell_full, type_pot = 'arith', scale = 1))
crit_larith <- c(list(point = crit_larith),
                  calc_cov_criterion(draws_logocc, draws_pot_larith, sig2 = sig2_larith,
                                     cell_full, type_pot = 'log_arith', scale = 1))

save(crit_arith, crit_larith, pred_occ, pred_pot_arith, pred_pot_larith, file = file.path(interim_results_dir,
                                                        'cv_total_density.Rda'))

