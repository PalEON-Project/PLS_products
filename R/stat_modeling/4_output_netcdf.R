## dimensions for x,y,draw_num
## variables for total, each taxon
## mean and sd to be derived from draws (check that mean is very similar to point fit from gam)

library(ncdf4)

output_netcdf_name <- paste0("PLS_biomass_western_v", product_version, ".nc")

x_grid <- sort(unique(pred_grid_west$x))
y_grid <- sort(unique(pred_grid_west$y))
x_res <- length(x_grid)
y_res <- length(y_grid)

make_albers_netcdf(name = 'biomass', longname = 'biomass, Mg per hectare', fn = output_netcdf_name, dir = output_dir, x = x_grid, y = y_grid, taxa = c('total', taxa), num_samples = n_stat_samples)

output_netcdf <- nc_open(file.path(output_dir, output_netcdf_name), write = TRUE)

mask <- nc_open(file.path(dataDir, 'paleonMask.nc'))
regions <- ncvar_get(mask, 'subregion', c(1,1),c(-1,-1))
wat <- ncvar_get(mask, 'water', c(1,1),c(-1,-1))
dom <- ncvar_get(mask, 'domain', c(1,1),c(-1,-1))

ncvar_put(output_netcdf, 'total', biomass_total$draws, start = c(1, 1, 1), count = c(x_res, y_res, num_samples))

for(p in seq_len(taxa)) 
    ncvar_put(output_netcdf, taxa[p], biomass_taxon[[taxa[p]]]$draws, start = c(1, 1, 1), count = c(x_res, y_res, num_samples))

nc_close(output_netcdf)
