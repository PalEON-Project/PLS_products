library(raster)
library(dplyr)
library(ncdf4)

if(shared_params_in_cell)
    load('point_with_biomass_shared.Rda') else load('point_with_biomass.Rda')

if(!'cell' %in% names(mw)) 
    mw <- mw %>% add_cells_to_dataset() 

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
