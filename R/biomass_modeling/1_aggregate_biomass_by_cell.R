##
# KH adding the code that has been used in the past to aggregate biomass, density, etc up to the grid cell:

spec.table <- read.csv(paste0('outputs/biomass_no_na_pointwise.ests_inilmi','_v',version, '.csv'))

#-------------------------Paleon gridding----------------------------------------------

# These are not the full tables since the include only the Paleon grid cells with points in the database.
#dcast rearranges the spec.table data by x, y and cell
count.table <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm=TRUE, value.var = 'count')

count.table.pt <- dcast(spec.table, Pointx + Pointy + cell ~ spec, sum, na.rm=TRUE, value.var = 'count')

unique.len <- function(x){length(unique(x))}

#melt data by by x, y, cell as rows and spec as columns and provide the sum or the number of unique points
biomass.trees  <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm=TRUE, value.var = 'count')

#biomass.trees represnts the number of trees in each category

biomass.points <- dcast(spec.table, x + y + cell ~ spec, unique.len, value.var = 'point')
#biomass.points represents the number of unique points
# This is to get the points with trees:
spec.adj <- spec.table
spec.adj$tree <- spec.table$spec %in% "No tree"

treed.points <- dcast(spec.adj, x + y + cell ~ tree, unique.len, value.var = 'point')
treed.points <- treed.points[order(treed.points$cell),]
colnames(treed.points) <- c('x', 'y', 'cell', 'tree', 'no tree')

# Now get the total number of plots per cell:
plots_per_cell <- dcast(spec.table, x + y + cell ~ ., unique.len, value.var = 'point')
plots_per_cell <- plots_per_cell[order(plots_per_cell$cell),]

points_per_cell_df <- data.frame(plots_per_cell,
                                 treed = treed.points$tree,
                                 untreed = treed.points$`no tree`)

colnames(points_per_cell_df)[4:5] <- c("total_pls_points", "total_treed_points")


# points by cell is the # of points 
points.by.cell <- rowSums(biomass.points[, 4:ncol(biomass.points)], na.rm=TRUE)
trees.by.cell  <- rowSums(biomass.trees[,!colnames(biomass.trees) %in% c('x', 'y', 'cell', 'No tree','Wet', 'Water')], na.rm=TRUE)

#calculate the sum of total density, basal area, biomass & diameter by cell and species
density.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'density')
basal.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'basal')
biomass.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'biom')
diam.table <-  dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'diams')



#calculate standard deviations of density, basal area, biomass, and diameters
density.sd.table <- dcast(spec.table, x + y  + cell ~ spec, sd, na.rm=TRUE, value.var = 'density')
basal.sd.table <- dcast(spec.table, x + y  + cell ~ spec, sd, na.rm=TRUE, value.var = 'basal')
diam.sd.table <-  dcast(spec.table, x + y  + cell ~ spec, sd, na.rm=TRUE, value.var = 'diams')
biom.sd.table <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm=TRUE, value.var = 'biom')

#  The function averages the estimates a single value for density, basal area, and biomass based on the points.by.cell
normalize <- function(x, mult = 2, value = points.by.cell) {
  x[,4:ncol(x)] <-  x[,4:ncol(x)] / value *mult; x}

density.table <- normalize(density.table)
basal.table <- normalize(basal.table)
diam.table <- normalize(diam.table, mult = 2.54, trees.by.cell)
biomass.table <- normalize(biomass.table)

# output preliminary density table to look at:
write.csv(density.table, "outputs/density.table_test.csv")

# get a density total value
density.table$total = rowSums(density.table[,4:ncol(density.table)], na.rm=TRUE)

mass.table <- merge(biomass.table, points_per_cell_df[,c("x", "y", "cell", "total_pls_points", "total_treed_points")], by = c("x", "y", "cell"))

#  We want rasterized versions of these tables with sums:
rast.fun <- function(x){
  
 to_grid <- data.frame(cell = x$cell, 
                  total = rowSums(x[,4:ncol(x)], na.rm=TRUE))
  
  empty <- rep(NA, ncell(base.rast))
  empty[to_grid$cell] <- to_grid$total
  setValues(base.rast, empty)
  }

base.rast <- raster(xmn = -71000, xmx = 2297000, ncols=296,
                    ymn = 58000,  ymx = 1498000, nrows = 180,
                    crs = '+init=epsg:3175')

#create total rasters for densi, basal area, diameters, biomass
dens     <- rast.fun(density.table)
count.up <- rast.fun(count.table)
basal    <- rast.fun(basal.table)
mdiam    <- rast.fun(diam.table); mdiam[mdiam==0] <- NA
biomass <- rast.fun(biomass.table)

writeRaster(biomass, "data/biomass.grd", overwrite=TRUE)
writeRaster(dens, "data/dens.grd", overwrite=TRUE)


#to get sd for these
rowset <- cbind(1:(nrow(spec.table)/2), (nrow(spec.table)/2+1):nrow(spec.table))

sd.table <- data.frame(cell = spec.table$cell[rowset[,1]],
                       density = spec.table$density[rowset[,1]] * 2,
                       basal   = spec.table$basal[rowset[,1]] + spec.table$basal[rowset[,2]],
                       biomass = spec.table$biom[rowset[,1]] + spec.table$biom[rowset[,2]],
                       diam = spec.table$diam[rowset[,1]] + spec.table$diam[rowset[,2]])

dens.sd <- dcast(sd.table, cell ~ ., fun.aggregate=sd, na.rm=TRUE, value.var = 'density')
basal.sd <- dcast(sd.table, cell ~ ., sd, na.rm=TRUE, value.var = 'basal')
biomass.sd <- dcast(sd.table, cell ~ ., sd, na.rm=TRUE, value.var = 'biomass')
diam.sd <- dcast(sd.table, cell ~., fun.aggregate = sd, na.rm = TRUE, value.var = 'diam')
colnames(dens.sd) <-c("cell", "dens.sd")
colnames(basal.sd) <-c("cell", "basal.sd")
colnames(biomass.sd) <-c("cell", "biomass.sd")
colnames(diam.sd) <-c("cell", "diam.sd")

#combine sd's with the tables
density.table <- merge(density.table, dens.sd, by = "cell")
diam.table <- merge(diam.table, diam.sd, by = "cell")



data.table <- data.frame(xyFromCell(dens, 1:ncell(dens)), 
                         stem.dens = getValues(dens),
                         basal.area = getValues(basal),
                         biomass = getValues(biomass),
                         diam = getValues(mdiam))

data.table <- data.table[!is.na(data.table[,3]), ] # remove places missing cells
write.csv(data.table, paste0('data/density.basal.biomass_alb','_v',1, '.csv'))


#  Now we need to add zero cells to the dataframe:
comp.table <- data.frame(xyFromCell(base.rast, 1:ncell(base.rast)), 
                         cell = 1:ncell(base.rast),
                         matrix(ncol = ncol(biomass.table)-3, nrow = ncell(base.rast)))

colnames(comp.table)[4:ncol(biomass.table)] <- colnames(biomass.table)[4:ncol(comp.table)]

reform <- function(x){
  comp.table[x$cell,] <- x
  comp.table
}


biomass.full <- reform(biomass.table)
density.full <- reform(density.table)
diameter.full <- reform(diam.table)
composition.table <- reform(basal.table)


#make composition table
composition.table <- composition.table[,!names(composition.table) %in% c('Water', "Wet", "No tree", "plss_pts")]
composition.table[,4:ncol(composition.table)] <- composition.table[,4:ncol(composition.table)]/rowSums(composition.table[,4:ncol(composition.table)], na.rm=TRUE)

# composition table without paleon grid NA values:
comp.inil <- composition.table[complete.cases(composition.table),]


#extract biomass.full table by the extent to get a datatable


biomass.df.test <- extract(biomass, extent(ind_il), cellnumbers=TRUE, xy=TRUE)
#biomass.df.test has the cell numbers of the extent of indiana, so I will use these

#get the biomass.full cells out of biomass.full
biomass.df.test.cells <- data.frame(biomass.df.test[,1])
colnames(biomass.df.test.cells) <- 'cell'
biomass.indiana <- merge(biomass.full, biomass.df.test.cells, by ="cell")

#to get biomass.pts and tree.pts, use merge(x, y, by, all=TRUE)
cells.xy.ind <- biomass.indiana[,1:3] #xy and cells in indiana
biomass.points.ind <- merge(cells.xy.ind, biomass.points, all=TRUE)
count.trees.ind <- merge(cells.xy.ind, count.table, all=TRUE)


add.v <- function(x, name, row.names){
  
  #  Quick file name formatter:
  
    p.ext <- '_alb'
  
  write.csv(x, paste0('data/outputs/', name,  '_v',version, '.csv'), row.names = row.names)
}

add.v(count.table, 'plss_counts', row.names=FALSE)
add.v(biomass.points, 'plss_points', row.names=FALSE)
add.v(count.table, "plss_composition", row.names=FALSE)
add.v(count.table, "plss_inil_composition", row.names=FALSE)
biomass.table$plsspts_cell <- points.by.cell

add.v(density.table, 'plss_density', row.names=FALSE)
add.v(basal.table, 'plss_basal', row.names=FALSE)
add.v(biomass.table, 'plss_biomass', row.names=FALSE)

add.v(biomass.full,  'plss_spec_biomass', row.names=FALSE) #full species biomass for indiana
add.v(density.full, 'plss_spec_density', row.names=FALSE)
add.v(diameter.full, 'plss_spec_diam', row.names=FALSE)
add.v(diam.table, 'plss_diam', row.names=FALSE)
