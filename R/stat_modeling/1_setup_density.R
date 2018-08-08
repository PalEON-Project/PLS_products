library(raster)
library(dplyr)
library(readr)
library(ncdf4)

load(file.path(interim_results_dir, 'point_with_density.Rda'))

if(!'cell' %in% names(mw)) 
    mw <- mw %>% add_paleon_grid()

taxa_conv <- read_csv(file.path(conversions_data_dir, level_3a_to_3s_conversion_file)) %>%
    rename(omit_western = "omit western")

taxa <- unique(taxa_conv$level3s[taxa_conv$omit_western == "no"])

## very few chestnut so lump them
taxa_conv <- taxa_conv %>% mutate(level3s = ifelse(level3s == 'Chestnut', "Other hardwood", level3s)) %>% select(level3a, level3s)

mw <- mw %>% 
    left_join(taxa_conv, by = c('L3_tree1' = 'level3a')) %>% rename(L3s_tree1 = level3s) %>%
    left_join(taxa_conv, by = c('L3_tree2' = 'level3a')) %>% rename(L3s_tree2 = level3s) %>%
    left_join(taxa_conv, by = c('L3_tree3' = 'level3a')) %>% rename(L3s_tree3 = level3s) %>%
    left_join(taxa_conv, by = c('L3_tree4' = 'level3a')) %>% rename(L3s_tree4 = level3s)



xy <- xyFromCell(base_raster, getValues(base_raster))
grid <- tibble(x = xy[,1], y = xy[,2], cell = getValues(base_raster))

mask <- nc_open(file.path(conversions_data_dir, 'paleonmask.nc'))
regions <- ncvar_get(mask, 'subregion', c(1,1),c(-1,-1))

x <- matrix(ncvar_get(mask, 'x', c(1),c(-1)), nrow = nrow(regions),
            ncol = ncol(regions))
y <- matrix(ncvar_get(mask, 'y', c(1),c(-1)),nrow = nrow(regions),
            ncol = ncol(regions), byrow = TRUE) 

west <- c(regions %in% paleon_regions_west)
paleon <- c(regions %in% paleon_regions)
pred_grid <- data.frame(x = c(x), y = c(y))
pred_grid_paleon <- pred_grid[paleon, ]
pred_grid_west <- pred_grid[west, ]
