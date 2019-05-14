library(raster)
library(dplyr)
library(readr)
library(ncdf4)

if(shared_params_in_cell) {
    load(file.path(interim_results_dir, 'point_with_biomass_shared.Rda'))
} else load(file.path(interim_results_dir, 'point_with_biomass.Rda'))

if(!'cell' %in% names(mw)) 
    mw <- mw %>% add_paleon_grid()

nr = nrow(mw)

## Can't calculate point-level biomass at points without density estimate.
mw <- mw %>% filter(!is.na(density_for_biomass))
cat("Excluding ", nr - nrow(mw), " points with missing density for biomass.\n")
## Can't calculate point-level biomass (though it would be low).
nr = nrow(mw)
mw <- mw %>% filter(!(num_trees == 1 & is.na(biomass1)))
cat("Excluding ", nr - nrow(mw), " one-tree points with missing biomass.\n")
## Could calculate for case with only one missing biomass, but it gets
## complicated to do in an unbiased way so exclude.
nr = nrow(mw)
mw <- mw %>% filter(!(num_trees == 2 & (is.na(biomass1) | is.na(biomass2))))
cat("Excluding ", nr - nrow(mw), " two-tree points with at least one missing biomass.\n")


taxa_conv <- read_csv(file.path(conversions_data_dir, level_3a_to_3s_conversion_file)) %>%
    rename(omit_western = "omit western")

taxa <- unique(taxa_conv$level3s[taxa_conv$omit_western == "no"])

## Note: in FIA work we do the taxa conversion at a different point; consider that here.

## Too few of certain taxa to treat separately
taxa_conv <- taxa_conv %>% mutate(level3s = ifelse(level3s %in% excluded_level3s, "Other hardwood", level3s)) %>% dplyr::select(level3a, level3s)

mw <- mw %>% 
    left_join(taxa_conv, by = c('L3_tree1' = 'level3a')) %>% rename(L3s_tree1 = level3s) %>%
    left_join(taxa_conv, by = c('L3_tree2' = 'level3a')) %>% rename(L3s_tree2 = level3s) %>%
    left_join(taxa_conv, by = c('L3_tree3' = 'level3a')) %>% rename(L3s_tree3 = level3s) %>%
    left_join(taxa_conv, by = c('L3_tree4' = 'level3a')) %>% rename(L3s_tree4 = level3s)
## Of course at this point tree3 and tree4 will never be used, so don't need to do last two lines.


## Overall PalEON grid with cell IDs.
xy <- xyFromCell(base_raster, getValues(base_raster))
grid <- tibble(x = xy[,1], y = xy[,2], cell = getValues(base_raster))

mask <- nc_open(file.path(conversions_data_dir, 'paleonmask.nc'))
regions <- ncvar_get(mask, 'subregion', c(1,1),c(-1,-1))

x <- matrix(ncvar_get(mask, 'x', c(1),c(-1)), nrow = nrow(regions),
            ncol = ncol(regions))
y <- matrix(ncvar_get(mask, 'y', c(1),c(-1)),nrow = nrow(regions),
            ncol = ncol(regions), byrow = TRUE) 

## full rectangular grid
pred_grid <- data.frame(x = c(x), y = c(y))
## subset to only land areas within PalEON states:
pred_grid_paleon <- pred_grid[regions %in% paleon_regions, ]
pred_grid_west <- pred_grid[regions %in% paleon_regions_west, ]

save(mw, taxa, pred_grid, pred_grid_paleon, pred_grid_west, grid, regions,
     file = file.path(interim_results_dir, 'cell_with_biomass_grid.Rda'))
