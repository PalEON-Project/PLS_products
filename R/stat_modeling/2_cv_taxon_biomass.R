## Fit statistical model to smooth the raw cell-level biomass via cross-validation
## to determine best upper-bound on amount of spatial smoothing.

## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

## This is very computationally-intensive and best done on a cluster.

library(dplyr)

warning("This code uses ~100 GB of RAM.")


load(file.path(interim_results_dir, paste0('cell_with_biomass_grid',
                                           ifelse(use_agb, '_agb', ''), '.Rda')))

if(use_mpi) {
    library(doMPI)
    cl <- startMPIcluster()
    registerDoMPI(cl)
} else {
    library(doParallel)
    if(n_cores == 0) {
        if(Sys.getenv("SLURM_JOB_ID") != "") {
            n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
        } else n_cores <- detectCores()
    }
    registerDoParallel(cores = n_cores)
}

taxa_to_fit <- taxa 
print(taxa_to_fit)

## set up CV folds
set.seed(1)
cells <- sample(unique(mw$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))
cell_full <- data.frame(cell = cells, fold = folds) %>% arrange(cell)

## Fit statistical model to each taxon and fold.
## Nested foreach will run separate tasks for each combination of taxon and fold.
output <- foreach(taxonIdx = seq_along(taxa_to_fit)) %:%
    foreach(i = seq_len(n_folds)) %dopar% {

        taxon <- taxa_to_fit[taxonIdx]
        ## add taxon-specific point-level biomass to dataset
        tmp <- mw %>% mutate(biomass_focal = calc_biomass_taxon(num_trees, biomass1, biomass2,
                                                                density_for_biomass, L3s_tree1, L3s_tree2, taxon))
        
        ## add total point-level biomass to dataset
        cell_full_taxon <- tmp %>% filter(!(is.na(biomass_focal))) %>% group_by(cell) %>% summarize(points_total = n())
        
        ## biomass stats averaged over occupied points
        cell_occ <- tmp %>% filter(!is.na(biomass_focal) & biomass_focal > 0) %>% 
            group_by(cell) %>%
            summarize(avg = mean(biomass_focal),
                      geom_avg = mean(log(biomass_focal)),
                      points_occ = n())
        
        ## should have total number of points, occupied number of points, and biomass stats (for occupied points)
        cell_full_taxon <- cell_full_taxon %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
            left_join(grid, by = c("cell" = "cell")) %>%
            mutate(points_occ = ifelse(is.na(points_occ), 0 , points_occ))
    
        ## get fold information
        cell_full_taxon <- cell_full %>% left_join(cell_full_taxon, by = c("cell" = "cell")) %>% arrange(cell)

        train <- cell_full_taxon %>% filter(fold != i)
        test <- cell_full_taxon %>% filter(fold == i) %>% arrange(cell)
        
        po <- fit(train, newdata = test, k_occ = k_occ_cv, unc = TRUE, use_bam = TRUE, bound_draws_low = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
        ppa1 <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', weight_scale = 1, unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
        ppl1 <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', weight_scale = 1, unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
        ppa70 <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', weight_scale = 70, unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
        ppl70 <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', weight_scale = 70, unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
        cat("taxon: ", taxonIdx, " ; fold: ", i, "\n", sep = "")
        list(po = po, ppa1 = ppa1, ppl1 = ppl1, ppa70 = ppa70, ppl70 = ppl70)
    }

pred_occ <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_occ_cv)))
pred_pot_arith1 <- pred_pot_larith1 <- sig2_arith1 <- sig2_larith1 <-
pred_pot_arith70 <- pred_pot_larith70 <- sig2_arith70 <- sig2_larith70 <-
    array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_pot_cv)))
    
dimnames(pred_occ)[[1]] <- taxa_to_fit
dimnames(pred_occ)[[3]] <- k_occ_cv
dimnames(pred_pot_arith1)[[3]] <- dimnames(pred_pot_larith1)[[3]] <- dimnames(sig2_arith1)[[3]] <- dimnames(sig2_larith1)[[3]] <-
dimnames(pred_pot_arith70)[[3]] <- dimnames(pred_pot_larith70)[[3]] <- dimnames(sig2_arith70)[[3]] <- dimnames(sig2_larith70)[[3]] <- k_pot_cv

draws_logocc <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_occ_cv), n_stat_samples))
draws_logpot_arith1 <- draws_logpot_larith1 <- draws_logpot_arith70 <- draws_logpot_larith70 <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_pot_cv), n_stat_samples))


for(taxonIdx in seq_along(taxa_to_fit))
    for(i in seq_len(n_folds)) {
        subn <- sum(cell_full$fold == i)
        pred_occ[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$po$pred_occ
        pred_pot_arith1[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$ppa1$pred_pot
        pred_pot_larith1[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$ppl1$pred_pot
        pred_pot_arith70[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$ppa70$pred_pot
        pred_pot_larith70[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$ppl70$pred_pot
        draws_logocc[taxonIdx, cell_full$fold == i, , ] <- output[[taxonIdx]][[i]]$po$draws_logocc
        draws_logpot_arith1[taxonIdx, cell_full$fold == i, , ] <- output[[taxonIdx]][[i]]$ppa1$draws_logpot
        draws_logpot_larith1[taxonIdx, cell_full$fold == i, , ] <- output[[taxonIdx]][[i]]$ppl1$draws_logpot
        draws_logpot_arith70[taxonIdx, cell_full$fold == i, , ] <- output[[taxonIdx]][[i]]$ppa70$draws_logpot
        draws_logpot_larith70[taxonIdx, cell_full$fold == i, , ] <- output[[taxonIdx]][[i]]$ppl70$draws_logpot
        sig2_arith1[taxonIdx, cell_full$fold == i, ] <- rep(output[[taxonIdx]][[i]]$ppa1$model_pot, each = subn)
        sig2_larith1[taxonIdx, cell_full$fold == i, ] <- rep(output[[taxonIdx]][[i]]$ppl1$model_pot, each = subn)
        sig2_arith70[taxonIdx, cell_full$fold == i, ] <- rep(output[[taxonIdx]][[i]]$ppa70$model_pot, each = subn)
        sig2_larith70[taxonIdx, cell_full$fold == i, ] <- rep(output[[taxonIdx]][[i]]$ppl70$model_pot, each = subn)
    }

## Assess results.

crit_point_arith1 <- crit_point_larith1 <- crit_point_arith70 <- crit_point_larith70 <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(crit_point_arith1)[[1]] <- dimnames(crit_point_larith1)[[1]] <- dimnames(crit_point_arith70)[[1]] <- dimnames(crit_point_larith70)[[1]] <- taxa_to_fit
dimnames(crit_point_arith1)[[2]] <- dimnames(crit_point_larith1)[[2]] <- dimnames(crit_point_arith70)[[2]] <- dimnames(crit_point_larith70)[[2]] <- k_occ_cv
dimnames(crit_point_arith1)[[3]] <- dimnames(crit_point_larith1)[[3]] <- dimnames(crit_point_arith70)[[3]] <- dimnames(crit_point_larith70)[[3]] <- k_pot_cv

crit_cov_arith1 <- crit_cov_larith1 <- crit_cov_arith70 <- crit_cov_larith70 <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(crit_cov_arith1)[[1]] <- dimnames(crit_cov_larith1)[[1]] <- dimnames(crit_cov_arith70)[[1]] <- dimnames(crit_cov_larith70)[[1]] <- taxa_to_fit
dimnames(crit_cov_arith1)[[2]] <- dimnames(crit_cov_larith1)[[2]] <- dimnames(crit_cov_arith70)[[2]] <- dimnames(crit_cov_larith70)[[2]] <- k_occ_cv
dimnames(crit_cov_arith1)[[3]] <- dimnames(crit_cov_larith1)[[3]] <- dimnames(crit_cov_arith70)[[3]] <- dimnames(crit_cov_larith70)[[3]] <- k_pot_cv

crit_length_arith1 <- crit_length_larith1 <- crit_length_arith70 <- crit_length_larith70 <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(crit_length_arith1)[[1]] <- dimnames(crit_length_larith1)[[1]] <- dimnames(crit_length_arith70)[[1]] <- dimnames(crit_length_larith70)[[1]] <- taxa_to_fit
dimnames(crit_length_arith1)[[2]] <- dimnames(crit_length_larith1)[[2]] <- dimnames(crit_length_arith70)[[2]] <- dimnames(crit_length_larith70)[[2]] <- k_occ_cv
dimnames(crit_length_arith1)[[3]] <- dimnames(crit_length_larith1)[[3]] <- dimnames(crit_length_arith70)[[3]] <- dimnames(crit_length_larith70)[[3]] <- k_pot_cv

crit_loglength_arith1 <- crit_loglength_larith1 <- crit_loglength_arith70 <- crit_loglength_larith70 <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(crit_loglength_arith1)[[1]] <- dimnames(crit_loglength_larith1)[[1]] <- dimnames(crit_loglength_arith70)[[1]] <- dimnames(crit_loglength_larith70)[[1]] <- taxa_to_fit
dimnames(crit_loglength_arith1)[[2]] <- dimnames(crit_loglength_larith1)[[2]] <- dimnames(crit_loglength_arith70)[[2]] <- dimnames(crit_loglength_larith70)[[2]] <- k_occ_cv
dimnames(crit_loglength_arith1)[[3]] <- dimnames(crit_loglength_larith1)[[3]] <- dimnames(crit_loglength_arith70)[[3]] <- dimnames(crit_loglength_larith70)[[3]] <- k_pot_cv

for(taxonIdx in seq_along(taxa_to_fit)) {
    ## extract raw data (again) for the taxon
    taxon <- taxa_to_fit[taxonIdx]
    ## add taxon-specific point-level biomass to dataset
    tmp <- mw %>% mutate(biomass_focal = calc_biomass_taxon(num_trees, biomass1, biomass2, density_for_biomass, L3s_tree1, L3s_tree2, taxon))
    
    ## add total point-level biomass to dataset
    cell_full_taxon <- tmp %>% filter(!(is.na(biomass_focal))) %>% group_by(cell) %>% summarize(points_total = n())
    
    ## biomass stats averaged over occupied points
    cell_occ <- tmp %>% filter(!is.na(biomass_focal) & biomass_focal > 0) %>% 
        group_by(cell) %>%
        summarize(avg = mean(biomass_focal),
                  geom_avg = mean(log(biomass_focal)),
                  points_occ = n())
    
    ## should have total number of points, occupied number of points, and biomass stats (for occupied points)
    cell_full_taxon <- cell_full_taxon %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
        left_join(grid, by = c("cell" = "cell")) %>%
        mutate(points_occ = ifelse(is.na(points_occ), 0 , points_occ)) %>% arrange(cell)

    cell_full_taxon <- cell_full %>% left_join(cell_full_taxon, by = c("cell" = "cell")) %>% arrange(cell)

    y <- cell_full_taxon$avg*cell_full_taxon$points_occ/cell_full_taxon$points_total  ## actual average biomass over all cells
    y[is.na(y)] <- 0  # cells with no points with trees (since $avg will be NA)
    y[is.na(cell_full_taxon$points_occ)] <- NA  # exclude cells with no valid points
    
    crit_point_arith1[taxonIdx, , ] <- calc_point_criterion(pred_occ[taxonIdx, , ], pred_pot_arith1[taxonIdx, , ],
                                                 cell_full_taxon$points_total, y, cv_max_biomass)
    crit_point_larith1[taxonIdx, , ] <- calc_point_criterion(pred_occ[taxonIdx, , ], pred_pot_larith1[taxonIdx, , ],
                                                       cell_full_taxon$points_total, y, cv_max_biomass)
    crit_point_arith70[taxonIdx, , ] <- calc_point_criterion(pred_occ[taxonIdx, , ], pred_pot_arith70[taxonIdx, , ],
                                                 cell_full_taxon$points_total, y, cv_max_biomass)
    crit_point_larith70[taxonIdx, , ] <- calc_point_criterion(pred_occ[taxonIdx, , ], pred_pot_larith70[taxonIdx, , ],
                                                       cell_full_taxon$points_total, y, cv_max_biomass)

    cell_full_taxon$obs <- y
    tmp <- calc_cov_criterion(draws_logocc[taxonIdx, , , ], draws_logpot_arith1[taxonIdx, , , ], sig2 = sig2_arith1[taxonIdx, , ],
                              cell_full_taxon, type_pot = 'arith', scale = 1)
    crit_cov_arith1[taxonIdx, , ] <- tmp$cov
    crit_length_arith1[taxonIdx, , ] <- tmp$length
    crit_loglength_arith1[taxonIdx, , ] <- tmp$loglength
    
    tmp <- calc_cov_criterion(draws_logocc[taxonIdx, , , ], draws_logpot_larith1[taxonIdx, , , ], sig2 = sig2_larith1[taxonIdx, , ],
                                     cell_full_taxon, type_pot = 'log_arith', scale = 1)
    crit_cov_larith1[taxonIdx, , ] <- tmp$cov
    crit_length_larith1[taxonIdx, , ] <- tmp$length
    crit_loglength_larith1[taxonIdx, , ] <- tmp$loglength

    tmp <- calc_cov_criterion(draws_logocc[taxonIdx, , , ], draws_logpot_arith70[taxonIdx, , , ], sig2 = sig2_arith70[taxonIdx, , ],
                                        cell_full_taxon, type_pot = 'arith', scale = 70)
    crit_cov_arith70[taxonIdx, , ] <- tmp$cov
    crit_length_arith70[taxonIdx, , ] <- tmp$length
    crit_loglength_arith70[taxonIdx, , ] <- tmp$loglength

    tmp <- calc_cov_criterion(draws_logocc[taxonIdx, , , ], draws_logpot_larith70[taxonIdx, , , ], sig2 = sig2_larith70[taxonIdx, , ],
                                     cell_full_taxon, type_pot = 'log_arith', scale = 70)
    crit_cov_larith70[taxonIdx, , ] <- tmp$cov
    crit_length_larith70[taxonIdx, , ] <- tmp$length
    crit_loglength_larith70[taxonIdx, , ] <- tmp$loglength

}

save(crit_point_arith1, crit_point_larith1, crit_point_arith70, crit_point_larith70,
     crit_cov_arith1, crit_cov_larith1, crit_cov_arith70, crit_cov_larith70,
     crit_length_arith1, crit_length_larith1, crit_length_arith70, crit_length_larith70,
     crit_loglength_arith1, crit_loglength_larith1, crit_loglength_arith70, crit_loglength_larith70,
     pred_occ, pred_pot_arith1, pred_pot_larith1, pred_pot_arith70, pred_pot_larith70,
     file = file.path(interim_results_dir, paste0('cv_taxon_biomass_k_low', ifelse(use_agb, '_agb',''), '.Rda')))


if(use_mpi) closeCluster(cl)
