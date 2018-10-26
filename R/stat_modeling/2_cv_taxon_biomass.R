## Fit statistical model to smooth the raw cell-level biomass via cross-validation
## to determine best upper-bound on amount of spatial smoothing.

## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

## This is very computationally-intensive and best done on a cluster.

library(dplyr)

load(file.path(interim_results_dir, 'cell_with_biomass_grid.Rda'))

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
        
        po <- fit(train, newdata = test, k_occ = k_occ_cv, unc = FALSE, use_bam = TRUE)
        ppa <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', unc = FALSE, use_bam = TRUE)
        ppl <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', unc = FALSE, use_bam = TRUE)
        cat("taxon: ", taxonIdx, " ; fold: ", i, "\n", sep = "")
        list(po$pred_occ, ppa$pred_pot, ppl$pred_pot)
    }

pred_occ <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_occ_cv)))
pred_pot_arith <- pred_pot_larith <- array(0, c(length(taxa_to_fit), nrow(cell_full), length(k_pot_cv)))
    
dimnames(pred_occ)[[1]] <- dimnames(pred_pot_arith)[[1]] <- dimnames(pred_pot_larith)[[1]] <- taxa_to_fit
dimnames(pred_occ)[[3]] <- k_occ_cv
dimnames(pred_pot_arith)[[3]] <- dimnames(pred_pot_larith)[[3]] <- k_pot_cv

for(taxonIdx in seq_along(taxa_to_fit))
    for(i in seq_len(n_folds)) {
        pred_occ[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]][[1]]
        pred_pot_arith[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]][[2]]
        pred_pot_larith[taxonIdx, cell_full$fold == i, ] <- output[[taxonIdx]][[i]][[3]]
    }

## Assess results.

critArith <- critLogArith <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(critArith)[[1]] <- dimnames(critLogArith)[[1]] <- taxa_to_fit
dimnames(critArith)[[2]] <- dimnames(critLogArith)[[2]] <- k_occ_cv
dimnames(critArith)[[3]] <- dimnames(critLogArith)[[3]] <- k_pot_cv

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
    
    critArith[taxonIdx, , ] <- calc_cv_criterion(pred_occ[taxonIdx, , ], pred_pot_arith[taxonIdx, , ],
                                                 cell_full_taxon$points_total, y, cv_max_biomass)
    critLogArith[taxonIdx, , ] <- calc_cv_criterion(pred_occ[taxonIdx, , ], pred_pot_larith[taxonIdx, , ],
                                                    cell_full_taxon$points_total, y, cv_max_biomass)
}

save(critArith, critLogArith, pred_occ, pred_pot_arith, pred_pot_larith,
     file = file.path(interim_results_dir, 'cv_taxon_biomass.Rda'))


if(use_mpi) closeCluster(cl)
