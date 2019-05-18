## Fit statistical model to smooth the raw cell level biomass.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

library(dplyr)
library(assertthat)

load(file.path(interim_results_dir, paste0('cell_with_biomass_grid',
                                           ifelse(use_agb, '_agb', ''), '.Rda')))

if(!exists('k_pot_total_biomass'))
    stop("Must specify 'k_pot_total_biomass'")
if(!exists('k_occ_total_biomass'))
    stop("Must specify 'k_occ_total_biomass'")

## note that in cases with one tree with 0 or NA diameter and therefore missing biomass
## we use the biomass for the other tree as the representative value
## this assumes that missingness of diameter is non-informative
cat("Found ", sum(mw$num_trees == 2 & (is.na(mw$biomass1) | is.na(mw$biomass2) |
                  mw$biomass1 == 0 | mw$biomass2 == 0)), " two-tree points with one tree missing biomass\n.")

biomass_avg <- mw %>% dplyr::select(biomass1, biomass2) %>% as.matrix(.) %>%
    apply(1, mean, na.rm = TRUE)

assert_that(sum(mw$num_trees==1 & is.na(mw$biomass1)) == 0,
            msg = "found 1-tree points with NA biomass")

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

## fit stats model
biomass_total <- fit(cell_full, newdata = pred_grid_west, k_occ = k_occ_total_biomass, k_pot = k_pot_total_biomass, return_model = TRUE, unc = TRUE, type_pot = fit_scale_biomass, num_draws = n_stat_samples, save_draws = TRUE, use_bam = TRUE)

save(biomass_total, file = file.path(interim_results_dir, paste0('fitted_total_biomass', ifelse(use_agb, '_agb',''), '.Rda')))
