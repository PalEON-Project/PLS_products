## Fit statistical model to smooth the raw cell level basal area.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average basal area for occupied points (called potential basal area).
## Estimated basal area is the product of occupancy and potential.

library(dplyr)

load(file.path(interim_results_dir, 'cell_with_basalarea_grid.Rda'))

## Allow for parallelization across taxa, including on Berkeley Statistics cluster with SLURM scheduler
library(doParallel)
if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}
registerDoParallel(cores = n_cores)

if(!exists('k_pot_taxon_biomass'))
    stop("Must specify 'k_pot_taxon_biomass'")
if(!exists('k_occ_taxon_biomass'))
    stop("Must specify 'k_occ_taxon_biomass'")

taxa_to_fit <- taxa 
print(taxa_to_fit)


basalarea_taxon <- foreach(taxonIdx = seq_along(taxa_to_fit)) %dopar% {
    taxon <- taxa_to_fit[taxonIdx]
    ## add taxon-specific point-level basal area to dataset
    tmp <- mw %>% mutate(basalarea_focal = calc_basalarea_taxon(num_trees, basalarea1, basalarea2, density_for_biomass, L3s_tree1, L3s_tree2, taxon))
    
    ## add total point-level basalarea to dataset
    cell_full_taxon <- tmp %>% filter(!(is.na(basalarea_focal))) %>% group_by(cell) %>% summarize(points_total = n())
    
    ## basal area stats averaged over occupied points
    cell_occ <- tmp %>% filter(!is.na(basalarea_focal) & basalarea_focal > 0) %>% 
        group_by(cell) %>%
    summarize(avg = mean(basalarea_focal),
              geom_avg = mean(log(basalarea_focal)),
              points_occ = n())
    
    ## should have total number of points, occupied number of points, and basal area stats (for occupied points)
    cell_full_taxon <- cell_full_taxon %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
        left_join(grid, by = c("cell" = "cell")) %>%
        mutate(points_occ = ifelse(is.na(points_occ), 0 , points_occ))
    
    ## fit stats model
    try(fit(cell_full_taxon, newdata = pred_grid_west, k_occ = k_occ_taxon_biomass, k_pot = k_pot_taxon_biomass, unc = TRUE, return_model = FALSE, type_pot = fit_scale_biomass, num_draws = n_stat_samples, save_draws = TRUE, use_bam = TRUE))
}

names(basalarea_taxon) <- taxa_to_fit
save(basalarea_taxon, file = file.path(interim_results_dir, 'fitted_taxon_basalarea.Rda'))
