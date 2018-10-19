## Fit statistical model to smooth the raw cell level biomass.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

library(dplyr)

load(file.path(interim_results_dir, 'cell_with_biomass_grid.Rda'))

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


biomass_taxon <- foreach(taxonIdx = seq_along(taxa_to_fit)) %dopar% {
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
        mutate(points_occ = ifelse(is.na(points_occ), 0 , points_occ))
    
    ## fit stats model
    try(fit(cell_full_taxon, newdata = pred_grid_west, k_occ = k_occ_taxon_biomass, k_pot = k_pot_taxon_biomass, unc = TRUE, return_model = FALSE, type_pot = fit_scale_biomass, num_draws = n_stat_samples, save_draws = TRUE, use_bam = TRUE))
}

names(biomass_taxon) <- taxa_to_fit
save(biomass_taxon, file = file.path(interim_results_dir, 'fitted_taxon_biomass.Rda'))
