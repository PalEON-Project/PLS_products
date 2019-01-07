## Outputs a netCDF file with dimensions of x, y, samples and
## variables for total basal area and basal area of each taxon

library(ncdf4)

load(file.path(interim_results_dir, paste0('fitted_total_basalarea.Rda')))
load(file.path(interim_results_dir, paste0('fitted_taxon_basalarea.Rda')))
load(file.path(interim_results_dir, 'cell_with_basalarea_grid.Rda')) # for grid info

taxa <- names(basalarea_taxon)
taxaNames <- gsub("/", ",", taxa)  ## netCDF interprets / as 'groups'

output_netcdf_name <- paste0("PLS_basalarea_western_v", product_version, ".nc")
output_netcdf_name_point <- paste0("PLS_basalarea_western_point_v", product_version, ".nc")

x_grid <- sort(unique(pred_grid$x))
y_grid <- sort(unique(pred_grid$y))

westernDomainX <- 1:146
westernDomainY <- 1:180

x_grid <- x_grid[westernDomainX]
y_grid <- y_grid[westernDomainY]
x_res <- length(x_grid)
y_res <- length(y_grid)

make_albers_netcdf(name = 'basalarea', longname = 'basal area, square meters per hectare,', units = 'm^2/ha', fn = output_netcdf_name, dir = output_dir, x = x_grid, y = y_grid, taxa = c('Total', taxaNames), num_samples = n_stat_samples)

make_albers_netcdf(name = 'basalarea', longname = 'basal area, square meters per hectare,', units = 'm^2/ha', fn = output_netcdf_name_point, dir = output_dir, x = x_grid, y = y_grid, taxa = c('Total', taxaNames), num_samples = 0)

output_netcdf <- nc_open(file.path(output_dir, output_netcdf_name), write = TRUE)
output_netcdf_point <- nc_open(file.path(output_dir, output_netcdf_name_point), write = TRUE)

mask <- nc_open(file.path(conversions_data_dir, 'paleonmask.nc'))
regions <- ncvar_get(mask, 'subregion', c(1,1),c(-1,-1))
wat <- ncvar_get(mask, 'water', c(1,1),c(-1,-1))
dom <- ncvar_get(mask, 'domain', c(1,1),c(-1,-1))

full <- array(as.numeric(NA), c(x_res, y_res, n_stat_samples))
full_point <- array(as.numeric(NA), c(x_res, y_res))

locs <- basalarea_total$locs
preds <- basalarea_total$pred[ , 1]
draws <- basalarea_total$draws

for(i in seq_len(nrow(locs))) {
    full_point[which(locs[i,1] == x_grid), which(locs[i,2] == y_grid)] <- preds[i]
    full[which(locs[i,1] == x_grid), which(locs[i,2] == y_grid), ] <- draws[i, ]
    if(i%%1000 == 0) print(i)
}

ncvar_put(output_netcdf, 'Total', full, start = c(1, 1, 1), count = c(x_res, y_res, n_stat_samples))
ncvar_put(output_netcdf_point, 'Total', full_point, start = c(1, 1), count = c(x_res, y_res))

for(p in seq_along(taxa)) {
    locs <- basalarea_taxon[[p]]$locs
    preds <- basalarea_taxon[[p]]$pred[ , 1]
    draws <- basalarea_taxon[[p]]$draws
    full <- array(as.numeric(NA), c(x_res, y_res, n_stat_samples))
    full_point <- array(as.numeric(NA), c(x_res, y_res))
    for(i in seq_len(nrow(locs))) {
        full_point[which(locs[i,1] == x_grid), which(locs[i,2] == y_grid)] <- preds[i]
        full[which(locs[i,1] == x_grid), which(locs[i,2] == y_grid), ] <- draws[i, ]
        if(i%%1000 == 0) print(i)
    }
    ncvar_put(output_netcdf, taxaNames[p], full, start = c(1, 1, 1), count = c(x_res, y_res, n_stat_samples))
    ncvar_put(output_netcdf_point, taxaNames[p], full_point, start = c(1, 1), count = c(x_res, y_res))
}

nc_close(output_netcdf)
nc_close(output_netcdf_point)



