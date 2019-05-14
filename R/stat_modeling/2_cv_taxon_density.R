## Fit statistical model to smooth the raw cell-level density via cross-validation
## to determine best upper-bound on amount of spatial smoothing.

## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.

library(dplyr)

load(file.path(interim_results_dir, 'cell_with_density_grid.Rda'))

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
        ## add taxon-specific point-level density to dataset
        tmp <- mw %>% mutate(density_focal = calc_density_taxon(num_trees, density, L3s_tree1, L3s_tree2, taxon))
        
        ## add total point-level density to dataset
        cell_full_taxon <- tmp %>% filter(!(is.na(density_focal))) %>% group_by(cell) %>% summarize(points_total = n())
        
        ## density stats averaged over occupied points
        cell_occ <- tmp %>% filter(!is.na(density_focal) & density_focal > 0) %>% 
            group_by(cell) %>%
            summarize(avg = mean(density_focal),
                      geom_avg = mean(log(density_focal)),
                      points_occ = n())
        
        ## should have total number of points, occupied number of points, and density stats (for occupied points)
        cell_full_taxon <- cell_full_taxon %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
            left_join(grid, by = c("cell" = "cell")) %>%
            mutate(points_occ = ifelse(is.na(points_occ), 0 , points_occ))
        
        ## get fold information
        cell_full_taxon <- cell_full %>% left_join(cell_full_taxon, by = c("cell" = "cell")) %>% arrange(cell)

        train <- cell_full_taxon %>% filter(fold != i)
        test <- cell_full_taxon %>% filter(fold == i)
        
        po <- fit(train, newdata = test, k_occ = k_occ_cv, unc = TRUE, use_bam = TRUE, bound_draws_low = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
        ppa <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
        ppl <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
        print(i, taxonIdx)
        list(po = po, ppa = ppa, ppl = ppl)
    }

pred_occ <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_occ_cv)))
pred_pot_arith <- pred_pot_larith <- sig2_arith <- sig2_larith <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_pot_cv)))
    
dimnames(pred_occ)[[1]] <- dimnames(pred_pot_arith)[[1]] <- dimnames(pred_pot_larith)[[1]] <- taxa_to_fit
dimnames(pred_occ)[[3]] <- k_occ_cv
dimnames(pred_pot_arith)[[3]] <- dimnames(pred_pot_larith)[[3]] <- dimnames(sig2_arith)[[3]] <-
    dimnames(sig2_larith)[[3]] <- k_pot_cv

draws_logocc <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_occ_cv), n_stat_samples))
draws_logpot_arith <- draws_logpot_larith <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_pot_cv), n_stat_samples))

for(taxonIdx in seq_along(taxa_to_fit))
    for(i in seq_len(n_folds)) {
        pred_occ[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$po$pred_occ
        pred_pot_arith[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$ppa$pred_pot
        pred_pot_larith[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$ppl$pred_pot
        draws_logocc[taxonIdx, cell_full$fold == i, , ] <- output[[taxonIdx]][[i]]$po$draws_logocc
        draws_logpot_arith[taxonIdx, cell_full$fold == i, , ] <- output[[taxonIdx]][[i]]$ppa$draws_logpot
        draws_logpot_larith[taxonIdx, cell_full$fold == i, , ] <- output[[taxonIdx]][[i]]$ppl$draws_logpot
        sig2_arith[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$ppa$model_pot
        sig2_larith[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]]$ppl$model_pot
    }

## Assess results.

crit_point_arith <- crit_point_larith <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(crit_point_arith)[[1]] <- dimnames(crit_point_larith)[[1]] <- taxa_to_fit
dimnames(crit_point_arith)[[2]] <- dimnames(crit_point_larith)[[2]] <- k_occ_cv
dimnames(crit_point_arith)[[3]] <- dimnames(crit_point_larith)[[3]] <- k_pot_cv

crit_cov_arith <- crit_cov_larith <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(crit_cov_arith)[[1]] <- dimnames(crit_cov_larith)[[1]] <- taxa_to_fit
dimnames(crit_cov_arith)[[2]] <- dimnames(crit_cov_larith)[[2]] <- k_occ_cv
dimnames(crit_cov_arith)[[3]] <- dimnames(crit_cov_larith)[[3]] <- k_pot_cv

crit_length_arith <- crit_length_larith <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(crit_length_arith)[[1]] <- dimnames(crit_length_larith)[[1]] <- taxa_to_fit
dimnames(crit_length_arith)[[2]] <- dimnames(crit_length_larith)[[2]] <- k_occ_cv
dimnames(crit_length_arith)[[3]] <- dimnames(crit_length_larith)[[3]] <- k_pot_cv

crit_loglength_arith <- crit_loglength_larith <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(crit_loglength_arith)[[1]] <- dimnames(crit_loglength_larith)[[1]] <- taxa_to_fit
dimnames(crit_loglength_arith)[[2]] <- dimnames(crit_loglength_larith)[[2]] <- k_occ_cv
dimnames(crit_loglength_arith)[[3]] <- dimnames(crit_loglength_larith)[[3]] <- k_pot_cv

for(taxonIdx in seq_along(taxa_to_fit)) {
    ## extract raw data (again) for the taxon
    taxon <- taxa_to_fit[taxonIdx]
    ## add taxon-specific point-level density to dataset
    tmp <- mw %>% mutate(density_focal = calc_density_taxon(num_trees, density, L3s_tree1, L3s_tree2, taxon))
    
    ## add total point-level density to dataset
    cell_full_taxon <- tmp %>% filter(!(is.na(density_focal))) %>% group_by(cell) %>% summarize(points_total = n())
    
    ## density stats averaged over occupied points
    cell_occ <- tmp %>% filter(!is.na(density_focal) & density_focal > 0) %>% 
        group_by(cell) %>%
        summarize(avg = mean(density_focal),
                  geom_avg = mean(log(density_focal)),
                  points_occ = n())
    
    ## should have total number of points, occupied number of points, and density stats (for occupied points)
    cell_full_taxon <- cell_full_taxon %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
        left_join(grid, by = c("cell" = "cell")) %>%
        mutate(points_occ = ifelse(is.na(points_occ), 0 , points_occ))
        
    y <- cell_full_taxon$avg*cell_full_taxon$points_occ/cell_full_taxon$points_total  ## actual average density over all cells
    y[is.na(y)] <- 0  # cells with no points with trees (since $avg will be NA)
    y[is.na(cell_full_taxon$points_occ)] <- NA  # exclude cells with no valid points (but for density, unlike biomass, this shouldn't occur because NA density values are excluded in 1_setup_density.R)

    crit_arith[taxonIdx, , ] <- calc_cv_criterion(pred_occ[taxonIdx, , ], pred_pot_arith[taxonIdx, , ], cell_full_taxon$points_total,
                                                 y, cv_max_density)
    crit_larith[taxonIdx, , ] <- calc_cv_criterion(pred_occ[taxonIdx, , ], pred_pot_larith[taxonIdx, , ], cell_full_taxon$points_total,
                                                   y, cv_max_density)

    cell_full_taxon$y <- y
    crit_arith <- c(list(point = crit_arith),
                     calc_cov_criterion(pred_occ[taxonIdx, , ], pred_pot_arith[taxonIdx, , ], sig2 = sig2_arith[taxonIdx, , ],
                                        cell_full_taxon, type_pot = 'arith', scale = 1))
    crit_larith <- c(list(point = crit_larith),
                      calc_cov_criterion(pred_occ[taxonIdx, , ], pred_pot_larith[taxonIdx, , ], sig2 = sig2_larith[taxonIdx, , ],
                                     cell_full_taxon, type_pot = 'log_arith', scale = 1))

}

save(crit_arith, crit_larith, pred_occ, pred_pot_arith, pred_pot_larith,
     file = file.path(interim_results_dir, 'cv_taxon_density.Rda'))


if(use_mpi) closeCluster(cl)


