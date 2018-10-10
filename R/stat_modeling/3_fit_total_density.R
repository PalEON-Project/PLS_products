## Fit statistical model to smooth the raw cell level density.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.

stop('need to determine k values and fit_scale for density and include info in config file')

load(file.path(interim_results_dir, 'cell_with_density_grid.Rda'))

if(!exists('k_pot_total_density'))
    stop("Must specify 'k_pot_total_density'")
if(!exists('k_occ_total_density'))
    stop("Must specify 'k_occ_total_density'")

## note that in cases with one tree with 0 or NA diameter and therefore missing density
## we use the density for the other tree as the representative value
## this assumes that missingness of diameter is non-informative
## 1384 2-tree points with tree2 diam missing
## 1317 2-tree points with tree1 diam missing

## total points per cell
cell_full <- mw %>% filter(!is.na(density)) %>%
    group_by(cell) %>% summarize(points_total = n())

## density stats averaged over occupied points
cell_occ <- mw %>% filter(!is.na(density) & density > 0) %>% 
    group_by(cell) %>%
    summarize(avg = mean(density),
              geom_avg = mean(log(density)),
              points_occ = n())

## note: given we fit stat model for potential density on log scale, can't really scale
## by count_occ unless we use geometric or arithmetic average, not log of arithmetic average
## but in cross-validation log of arithmetic average seems to work better for biomass
## so assuming that is the case for density for the moment as well

## should have 'total', 'occupied', and density stats (for occupied points)
cell_full <- cell_full %>% left_join(cell_occ, by = c("cell" = "cell")) %>%
    left_join(grid, by = c("cell" = "cell")) %>%
    mutate(points_occ = ifelse(is.na(points_occ), 0, points_occ))

## fit stats model
density_total <- fit(cell_full, newdata = pred_grid_west, k_occ = k_occ_total_density, k_pot = k_pot_total_density, return_model = TRUE, unc = TRUE, type_pot = fit_scale_density, num_draws = n_stat_samples, save_draws = TRUE, use_bam = TRUE)

save(density_total, file = file.path(interim_results_dir, 'fitted_total_density.Rda'))
