library(raster)
library(dplyr)
library(ncdf4)

if(!exists('k_pot_total'))
    stop("Must specify 'k_pot_total'")
if(!exists('k_occ_total'))
    stop("Must specify 'k_occ_total'")

biomass_avg <- mw %>% dplyr::select(biomass1, biomass2) %>% as.matrix(.) %>%
    apply(1, mean, na.rm = TRUE)

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
## by count_occ unless we use geometric avg

## should have 'total', 'occupied', and biomass stats (for occupied points)
cell_full <- cell_full %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
    left_join(grid, by = c("cell" = "cell")) %>%
    mutate(points_occ = ifelse(is.na(points_occ), 0, points_occ))

## fit stats model
biomass_total <- fit(cell_full, newdata = pred_grid_west, k_occ = k_occ_total, k_pot = k_pot_total,
                     unc = TRUE, type_pot = 'log_arith', num_draws = 1000)
