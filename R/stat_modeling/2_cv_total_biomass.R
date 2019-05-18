## Fit statistical model to smooth the raw cell-level biomass via cross-validation
## to determine best upper-bound on amount of spatial smoothing.

## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

## Note that in cases with one tree with 0 or NA diameter and therefore missing biomass
## we use the biomass for the other tree as the representative value
## this assumes that missingness of diameter is non-informative

library(dplyr)

load(file.path(interim_results_dir, paste0('cell_with_biomass_grid',
                                           ifelse(use_agb, "_agb", ""), ".Rda")))

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
pred_pot_arith1 <- pred_pot_larith1 <- pred_pot_arith70 <- pred_pot_larith70 <-
    sig2_arith1 <- sig2_arith70 <- sig2_larith1 <- sig2_larith70 <- matrix(0, nrow(cell_full), length(k_pot_cv))
dimnames(pred_pot_arith1)[[2]] <- dimnames(pred_pot_larith1)[[2]] <-
    dimnames(pred_pot_arith70)[[2]] <- dimnames(pred_pot_larith70)[[2]] <-
    dimnames(sig2_arith1)[[2]] <- dimnames(sig2_arith70)[[2]] <-
    dimnames(sig2_larith1)[[2]] <- dimnames(sig2_larith70)[[2]] <- k_pot_cv

draws_logocc <- array(0, c(nrow(cell_full), length(k_occ_cv), n_stat_samples))
draws_logpot_arith1 <- draws_logpot_larith1 <- draws_logpot_arith70 <- draws_logpot_larith70 <- array(0, c(nrow(cell_full), length(k_pot_cv), n_stat_samples))

n_folds <- max(cell_full$fold)
output <- foreach(i = seq_len(n_folds)) %dopar% {
    train <- cell_full %>% filter(fold != i)
    test <- cell_full %>% filter(fold == i)
    
    po <- fit(train, newdata = test, k_occ = k_occ_cv, unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
    ppa1 <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', weight_scale = 1, unc = TRUE,
                use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
    ppa70 <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', weight_scale = 70, unc = TRUE,
                 use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
    ppl1 <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', weight_scale = 1, unc = TRUE,
                use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
    ppl70 <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', weight_scale = 70, unc = TRUE,
                 use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
    cat("n_fold: ", i, " ", date(), "\n")
    list(po = po, ppa1 = ppa1, ppa70 = ppa70, ppl1 = ppl1, ppl70 = ppl70)
}


for(i in seq_len(n_folds)) {
    subn <- sum(cell_full$fold == i)
    pred_occ[cell_full$fold == i, ] <- output[[i]]$po$pred_occ
    pred_pot_arith1[cell_full$fold == i, ] <- output[[i]]$ppa1$pred_pot
    pred_pot_larith1[cell_full$fold == i, ] <- output[[i]]$ppl1$pred_pot
    pred_pot_arith70[cell_full$fold == i, ] <- output[[i]]$ppa70$pred_pot
    pred_pot_larith70[cell_full$fold == i, ] <- output[[i]]$ppl70$pred_pot

    draws_logocc[cell_full$fold == i, , ] <- output[[i]]$po$draws_logocc
    draws_logpot_arith1[cell_full$fold == i, , ] <- output[[i]]$ppa1$draws_logpot
    draws_logpot_larith1[cell_full$fold == i, , ] <- output[[i]]$ppl1$draws_logpot
    draws_logpot_arith70[cell_full$fold == i, , ] <- output[[i]]$ppa70$draws_logpot
    draws_logpot_larith70[cell_full$fold == i, , ] <- output[[i]]$ppl70$draws_logpot
    ## when model not returned, we do get back the 'sig2' value as the model_pot
    sig2_arith1[cell_full$fold == i, ] <- rep(output[[i]]$ppa1$model_pot, each = subn)
    sig2_larith1[cell_full$fold == i, ] <- rep(output[[i]]$ppl1$model_pot, each = subn)
    sig2_arith70[cell_full$fold == i, ] <- rep(output[[i]]$ppa70$model_pot, each = subn)
    sig2_larith70[cell_full$fold == i, ] <- rep(output[[i]]$ppl70$model_pot, each = subn)
}



## Assess results

y <- cell_full$avg*cell_full$points_occ/cell_full$points_total ## actual average biomass over all cells
y[is.na(y)] <- 0   # cells with no points with trees (since $avg will be NA)

crit_arith1 <- calc_point_criterion(pred_occ, pred_pot_arith1, cell_full$points_total,
                               y, cv_max_biomass)
crit_larith1 <- calc_point_criterion(pred_occ, pred_pot_larith1, cell_full$points_total,
                                  y, cv_max_biomass)
crit_arith70 <- calc_point_criterion(pred_occ, pred_pot_arith70, cell_full$points_total,
                               y, cv_max_biomass)
crit_larith70 <- calc_point_criterion(pred_occ, pred_pot_larith70, cell_full$points_total,
                                  y, cv_max_biomass)

cell_full$obs <- y
crit_arith1 <- c(list(point = crit_arith1),
                 calc_cov_criterion(draws_logocc, draws_logpot_arith1, sig2 = sig2_arith1,
                                    cell_full, type_pot = 'arith', scale = 1))
crit_larith1 <- c(list(point = crit_larith1),
                  calc_cov_criterion(draws_logocc, draws_logpot_larith1, sig2 = sig2_larith1,
                                     cell_full, type_pot = 'log_arith', scale = 1))
crit_arith70 <- c(list(point = crit_arith70),
                  calc_cov_criterion(draws_logocc, draws_logpot_arith70, sig2 = sig2_arith70,
                                     cell_full, type_pot = 'arith', scale = 70))
crit_larith70 <- c(list(point = crit_larith70),
                   calc_cov_criterion(draws_logocc, draws_logpot_larith70, sig2 = sig2_larith70,
                                      cell_full, type_pot = 'log_arith', scale = 70))

save(crit_arith1, crit_larith1, crit_arith70, crit_larith70,
     pred_occ, pred_pot_arith1, pred_pot_larith1, pred_pot_arith70, pred_pot_larith70, file = file.path(interim_results_dir,
                                                        paste0('cv_total_biomass', ifelse(use_agb, '_agb',''), '.Rda')))
