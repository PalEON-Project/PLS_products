## Fit statistical model to smooth the raw cell level density.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.

library(dplyr)

load(file.path(interim_results_dir, 'cell_with_density_grid.Rda'))

## Allow for parallelization across taxa, including on Berkeley Statistics cluster with SLURM scheduler
library(doParallel)
if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}
registerDoParallel(cores = n_cores)

if(!exists('k_pot_taxon_density'))
    stop("Must specify 'k_pot_taxon_density'")
if(!exists('k_occ_taxon_density'))
    stop("Must specify 'k_occ_taxon_density'")

taxa_to_fit <- taxa 
print(taxa_to_fit)

density_taxon <- foreach(taxonIdx = seq_along(taxa_to_fit)) %dopar% {
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
    
    ## fit stats model
    try(fit(cell_full_taxon, newdata = pred_grid_west, k_occ = k_occ_taxon_density, k_pot = k_pot_taxon_density, unc = TRUE, return_model = FALSE, type_pot = fit_scale_density, num_draws = n_stat_samples, save_draws = TRUE, use_bam = TRUE, bound_draws_low = TRUE))
}

names(density_taxon) <- taxa_to_fit
save(density_taxon, file = file.path(interim_results_dir, 'fitted_taxon_density.Rda'))
