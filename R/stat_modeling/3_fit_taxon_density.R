## Fit statistical model to smooth the raw point level density.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.

## Run-time (without cross-validation) approximately XYZ with k values of occ: 2500, pot: 2000

if(!exists('k_pot_taxon'))
    stop("Must specify 'k_pot_taxon'")
if(!exists('k_occ_taxon'))
    stop("Must specify 'k_occ_taxon'")

taxa_to_fit <- taxa # or put a single taxon of interest here

print(taxa_to_fit)

if(Sys.getenv("SLURM_JOB_ID") != "")
    nCores <- Sys.getenv("SLURM_NTASKS")

library(doParallel)
registerDoParallel(cores = nCores)
density_taxon <- foreach(taxonIdx = seq_along(taxa_to_fit)) %dopar% {
    taxon <- taxa[taxonIdx]
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
    
    if(do_cv) {
        
        k_occ <- c(100,250,500,1000,1500,2000,2500)
        k_pot = c(100,250,500,1000,1500,2000,2500,3000,3500)
        
        set.seed(1)
        cells <- sample(unique(cell_full_taxon$cell), replace = FALSE)
        folds <- rep(1:n_folds, length.out = length(cells))
        
        cell_full_taxon <- cell_full_taxon %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

        stop("this cv code won't work yet")
        density_taxon <- fit_cv(cell_full_taxon, k_occ, k_pot, n_cores)
        
        ## assess results
        y <- cell_full_taxon$avg*cell_full_taxon$count/cell_full_taxon$total  ## actual average density over all cells
        y[is.na(y)] <- 0
        y[y > mx] <- mx
        
        
        critArith <- calc_cv_criterion(density_taxon$pred_occ, density_taxon$pred_pot_arith, cell_full_taxon$points_total,
                                       cell_full_taxon$avg*cell_full_taxon$count/cell_full_taxon$total, 200)
        critLogArith <- calc_cv_criterion(density_taxon$pred_occ, density_taxon$pred_pot_larith, cell_full_taxon$points_total,
                                          cell_full_taxon$avg*cell_full_taxon$count/cell_full_taxon$total, 200)

    } else {
        ## fit stats model
        output <- try(fit(cell_full_taxon, newdata = pred_grid_west, k_occ = k_occ_taxon, k_pot = k_pot_taxon, unc = TRUE, return_model = FALSE, type_pot = 'log_arith', num_draws = n_stat_samples, save_draws = TRUE))
        output$critArith <- output$critLogArith <- NULL
        ## save(density_taxon, file = file.path(interim_results_dir, 'fitted_taxon_density2.Rda'))
        ## cat("Finished taxon: ", taxon, "\n")
    }
    output
}

names(density_taxon) <- taxa_to_fit
save(density_taxon, file = file.path(interim_results_dir, 'fitted_taxon_density.Rda'))
