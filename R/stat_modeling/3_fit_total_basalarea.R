## Fit statistical model to smooth the raw cell level basal area.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average basal area for occupied points (called potential basal area).
## Estimated basal area is the product of occupancy and potential.

library(dplyr)

load(file.path(interim_results_dir, 'cell_with_basalarea_grid.Rda'))

if(!exists('k_pot_total_biomass'))
    stop("Must specify 'k_pot_total_biomass'")
if(!exists('k_occ_total_biomass'))
    stop("Must specify 'k_occ_total_biomass'")

## note that in cases with one tree with 0 or NA diameter and therefore missing basal area
## we use the basal area for the other tree as the representative value
## this assumes that missingness of diameter is non-informative
## 1384 2-tree points with tree2 diam missing
## 1317 2-tree points with tree1 diam missing

basalarea_avg <- mw %>% dplyr::select(basalarea1, basalarea2) %>% as.matrix(.) %>%
    apply(1, mean, na.rm = TRUE)

## add total point-level basalarea to dataset
mw <- mw %>% mutate(basalarea = ifelse(num_trees == 0, 0,
                              ifelse(num_trees == 1, basalarea1, basalarea_avg)))

## total points per cell
cell_full <- mw %>% filter(!(is.na(basalarea) | is.na(density_for_biomass))) %>%
    group_by(cell) %>% summarize(points_total = n())

## basalarea stats averaged over occupied points
cell_occ <- mw %>% filter(!(is.na(basalarea) | is.na(density_for_biomass)) & basalarea > 0) %>% 
    group_by(cell) %>%
    summarize(avg = mean(basalarea*density_for_biomass),
              geom_avg = mean(log(basalarea*density_for_biomass)),
              points_occ = n())

## note: given we fit stat model for potential basal area on log scale, can't really scale
## by count_occ unless we use geometric or arithmetic average, not log of arithmetic average
## but in cross-validation log of arithmetic average seems to work better

## should have 'total', 'occupied', and basalarea stats (for occupied points)
cell_full <- cell_full %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
    left_join(grid, by = c("cell" = "cell")) %>%
    mutate(points_occ = ifelse(is.na(points_occ), 0, points_occ))

## fit stats model
basalarea_total <- fit(cell_full, newdata = pred_grid_west, k_occ = k_occ_total_biomass, k_pot = k_pot_total_biomass, return_model = TRUE, unc = TRUE, type_pot = fit_scale_biomass, num_draws = n_stat_samples, save_draws = TRUE, use_bam = TRUE)

save(basalarea_total, file = file.path(interim_results_dir, 'fitted_total_basalarea.Rda'))
