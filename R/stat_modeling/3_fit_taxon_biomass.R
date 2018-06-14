## Fit statistical model to smooth the raw point level biomass.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

## Run-time (without cross-validation) approximately with k values of occ: 2500, pot: 2000

if(!exists('k_pot_taxon'))
    stop("Must specify 'k_pot_taxon'")
if(!exists('k_occ_taxon'))
    stop("Must specify 'k_occ_taxon'")

taxa_to_fit <- taxa # or put a single taxon of interest here

biomass_taxon <- critArith <- critLogArith <- list()
length(biomass_taxon) <- length(critArith) <- length(critLogArith) <- length(taxa_to_fit)

for(taxon in taxa_to_fit) {
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
    
    if(do_cv) {
        
        k_occ <- c(100,250,500,1000,1500,2000,2500)
        k_pot = c(100,250,500,1000,1500,2000,2500,3000,3500)
        
        set.seed(1)
        cells <- sample(unique(cell_full_taxon$cell), replace = FALSE)
        folds <- rep(1:n_folds, length.out = length(cells))
        
        cell_full_taxon <- cell_full_taxon %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))
        
        results <- fit_cv(cell_full_taxon, k_occ, k_pot, n_cores)
        
        ## assess results
        y <- cell_full_taxon$avg*cell_full_taxon$count/cell_full_taxon$total  ## actual average biomass over all cells
        y[is.na(y)] <- 0
        y[y > mx] <- mx
        
        
        critArith[[taxon]] <- calc_cv_criterion(results$pred_occ, results$pred_pot_arith, cell_full_taxon$points_total,
                                       cell_full_taxon$avg*cell_full_taxon$count/cell_full_taxon$total, 200)
        critLogArith[[taxon]] <- calc_cv_criterion(results$pred_occ, results$pred_pot_larith, cell_full_taxon$points_total,
                                          cell_full_taxon$avg*cell_full_taxon$count/cell_full_taxon$total, 200)
        
    }
    
    
    ## fit stats model
    biomass_taxon[[taxon]] <- try(fit(cell_full_taxon, newdata = pred_grid_west, k_occ = k_occ_taxon, k_pot = k_pot_taxon, unc = TRUE, return_model = TRUE, type_pot = 'log_arith', num_draws = n_stat_samples))
    save(biomass_taxon, file = file.path(interim_results_dir, 'fitted_taxon_biomass.Rda'))
    cat("Finished taxon: ", taxon, "\n")
}


