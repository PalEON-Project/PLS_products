## to be based on KAH 01b_calculate_full_density.r (and 02_calculate_density.R?)
## calculate point-level density estimates

library(readr)
library(dplyr)


load('cleaned_point.Rda')

## determine number of trees per point

num_trees <- rep(2, nrow(mw))
## only 1 NA in L3_tree1; check back on where it came from
num_trees[is.na(mw$L3_tree1) | mw$L3_tree1 == "No tree"] <- 0
num_trees[!is.na(mw$L3_tree1) & mw$L3_tree1 != "No tree" & (mw$L3_tree2 == "No tree" | is.na(mw$L3_tree2))] <- 1

## 1714 locations in MN with one tree and water for 2nd tree; presumably treat this as
## one tree given other locations would be categorized as 1-tree
## note these locations overlap pretty well with spatial distribution of other 1-tree points
num_trees[mw$L3_tree1 != "No tree" & mw$L3_tree2 == "Water"] <- 1

## ~600 points in northern MN with no info that are classified as 0-tree
## the '_' vegtype in MN has 900 points many in northern MN and including some straight
## lines, suggesting we may want to throw these out

mw <- mw %>% mutate(num_trees = num_trees)

## Apply various exclusion criteria

## remove 1-tree 0-distance (or unknown distance) points as unclear what to use for density
## there are clusters of such points in LP of MI, southern IL and southern IN, and around Green Bay
## 3413 points
mw <- mw %>% filter(!(num_trees == 1 & (is.na(dist1) | dist1 == 0)))

## can't calculate density for points with missing distances: 1683 points
mw <- mw %>% filter(!(num_trees == 2 & (is.na(dist1) | is.na(dist2))))
## can't calculate biomass for points with missing diameters: 3338 points
## not doing this for the moment as these trees might be omitted when we also omit trees < 8 inches
## mw <- mw %>% filter(!(num_trees == 2 & (is.na(diam1) | is.na(diam2) | diam1 == 0 | diam2 == 0)))

if(any(mw[ , paste0('az', 1:4)] == 360)) {
    warning("Found some azimuths = 360")
    mw <- mw %>% mutate(az1 = ifelse(az1 == 360, 0, az1),
                        az2 = ifelse(az2 == 360, 0, az2),
                        az3 = ifelse(az3 == 360, 0, az3),
                        az4 = ifelse(az4 == 360, 0, az4))
}


## about 50 trees > 100 in diameter: 42 points
mw <- mw %>% filter(!( (!is.na(diam1) & diam1 > 100) | (!is.na(diam2) & diam2 > 100) ))

## note allometries mostly wouldn't go above 80 cm = 32 inches

## Plans for tricky cases:

## keep 2-tree points regardless of distances and truncate density (at say 1000 for now)
## use 2-tree points with small or NA diameter trees for density calculation,
## but treat trees with missing diam or diam < 8 as missing in terms of biomass calcs
## when both trees at a point are missing, that is simple: omit the point since
## the surveyor did not provide the info we need - two trees >= 8 are presumably there
## but we don't know anything about them
## when one tree at a point is missing, it is tricky because the missing tree could be
## of any taxon but we don't want to impute a zero for everything other than the one tree
## as that is biased low overall since we know the missing tree is some taxon
## I think we may need to omit these points too

## keep 1-tree points and set density to something like 1 tree per unit circle of radius (150 links for now)
## ~3.5 trees/ha
## treat trees with missing diam or small diameter as biomass of zero since we know biomass is low

if(F) {  # temporary exploratory code

mw1 = mw[mw$num_trees == 1,]
mw2 = mw[mw$num_trees == 2,]

## small tree points
small <- mw2 %>% filter(mw2$diam1 < 8 | mw2$diam2 < 8)

if(any(is.na(mw$L3_tree1[mw$num_trees == 1])))
   stop("Missing taxon for some 1-tree points.")

mw2$az1[mw2$az1 == 360] = 0
mw2$az2[mw2$az2 == 360] = 0
q1=floor(mw2$az1/90)
q2=floor(mw2$az2/90)
}

# about 10k points with points known to be in same azimuth; about 11k with missing az1, 8500 missing az2

##  - deal with azimuth for zero-distance trees - probably fine to have either numeric or NA values
## as will simply upperbound these densities


## Sensitivity analyses:
## w/ and w/o 1-tree points - some increase w/o 1-tree points but seems not generally more than 15-20% per cell
## perhaps do analysis where use the <8 inch trees to biomass estimates?


## Correction factors to account for surveyor sampling 'design'
corr_factors <- read_csv(file.path(conversions_data_dir, correction_factors_file))

## fix state names to match tree data
names_df <- data.frame(state = c('IN','IL','Michigan','Minnesota','Wisconsin'),
                       new_state = c('IN','IL','MI','MN','WI'), stringsAsFactors = FALSE)
corr_factors <- corr_factors %>% left_join(names_df, by = c('state' = 'state')) %>%
    dplyr::select(-state) %>% rename(state = new_state)

## TODO: remove when pull in updated corrections file
corr_factors <- corr_factors %>% mutate(sectioncorner = ifelse(sectioncorner == 'quartersection', 'quarter-section', sectioncorner))

## TODO: confirm with Charlie that 2nQ factors are never used
mw <- mw %>% mutate(density = calc_stem_density(mw, corr_factors))

save(mw, file = 'point_with_density.Rda')


if(FALSE) {  # Simon/Kelly older code

# calculates basal area and stem density using morisita function:
estimates <- morisita(final.data, correction.factor, veil = TRUE)

stem.density <- estimates[[1]]
basal.area <- estimates[[2]]

# there are some very high estimates of stem density 
summary(stem.density)
summary(basal.area)
zero.trees <- is.na(stem.density) 

#set stem.density where there are zero trees due to No tree or Wet or Water to 0
stem.density[zero.trees] <- 0
basal.area[zero.trees] <- 0

summary(stem.density)
summary(basal.area)

# make into a data frame and export as csv
stem.density <- data.frame(stem.density, basal.area, final.data)
write.csv(stem.density, paste0('outputs/IN_IL_MI_densestimates_v',version,'.csv'))



#----------------------------Density Regridding------------------------------

##need to regrid the density estimates onto the paleon centroids
##create base raster that is extent of midwest domain

base.rast <- raster(xmn = -71000, xmx = 2297000, ncols=296,
                  ymn = 58000,  ymx = 1498000, nrows = 180,
                    crs = '+init=epsg:3175')



#coordinates(final.data)<- ~PointX+PointY

#create spatial object with density, basal area & diameters data
stem.density <- data.frame(x = final.data$PointX, 
                           y = final.data$PointY,
                           corner = final.data$corner,
                           density = stem.density$stem.density,
                           basal   = stem.density$basal.area,
                           state = final.data$state,
                           township = final.data$Township)#,
                          #diams = rowMeans(diams[,1:2], na.rm=TRUE) * 2.54)


# ---------------------fixing some lingering data naming issues:-------------------


#fix the captalized "No tree" problem
species[species == 'No Tree'] <- 'No tree'
species[species==""]<- "No tree"

#change all No tree densities to 0
stem.density$density[species[,1] == 'No tree'| species[,2]=='No tree'] <- 0
#classify trees as zero or as wet trees
zero.trees <- is.na(stem.density$density) & (species[,2] %in% c('No tree') | species[,1] %in% c('No tree'))

#designate all zero trees as density of 0
stem.density$density[zero.trees] <- 0
stem.density$basal[zero.trees] <- 0


# kill cells with na for x or y:
stem.density <- stem.density[!is.na(stem.density$x),]

# make stem.density spatial
coordinates(stem.density) <- ~x+y
proj4string(stem.density)<-CRS('+init=epsg:3175')

# write to an arcGIS compatible shapefile
writeOGR(obj = stem.density, dsn = "outputs/stem_density_alb_v1.7-5.shp", layer = "stem_density_alb_v1.7-5", driver = "ESRI Shapefile", overwrite=TRUE)


#------------------------Formatting for biomass estimation-------------------------


numbered.rast <- setValues(base.rast, 1:ncell(base.rast))
numbered.cell <- raster::extract(numbered.rast, spTransform(stem.density,CRSobj=CRS('+init=epsg:3175')))

final.data <- data.frame(final.data)
#final.data <- read.csv(paste0("outputs/ndilinpls_for_density_v",version,".csv"), stringsAsFactors = FALSE)

#create dataframe with stem density, speceies
spec.table <- data.frame(PointX = final.data$PointX, 
                         PointY = final.data$PointY,
                         cell = numbered.cell,
                         spec = c(as.character(final.data$species1),as.character(final.data$species2)),
                         count = 1,
                         point = 1:nrow(final.data),
                         density = rep(stem.density$density/2, 2),
                         basal =  rep(stem.density$basal/2, 2),
                         diams = c(final.data$diam1, final.data$diam2),
                         dists = c(final.data$dist1, final.data$dist2),
                         state = final.data$state,
                         corner = final.data$corner, 
                         township = final.data$Township,
                         stringsAsFactors = FALSE)



                 
# This section should be moved to estimating biomss
#-----------------Estimating Biomass from density and diameter-------------------
# changing column names
spec.table$Pointx <- spec.table$PointX
spec.table$Pointy <- spec.table$PointY
spec.table[,1:2] <- xyFromCell(base.rast, spec.table$cell)
colnames(spec.table)[1:2] <- c("x", "y")
# read in table with allometric equations for each taxa
biom.table <- read.csv('data/plss.pft.conversion_v0.1-1.csv', 
                       stringsAsFactors = FALSE)

# this function calculates biomass of an individual tree using taxa-specific allometric equations
form <- function(x) {
  
  eqn <- match(x[c('spec')], biom.table[,1])
  eqn[is.na(eqn)] <- 1  #  Sets it up for non-tree.
  
  b0 <- biom.table[eqn,2]
  b1 <- biom.table[eqn,3]
  
  biomass <- exp(b0 + b1 * log(as.numeric(x[c('diams')])))
  biomass
  
}

#  This is the biomass of individual trees.  It needs to be converted into
#  a stand level value, through the stem density estimate  The values are
#  in kg.

#biomass <- rep(NA, nrow(spec.table))
biomass <- apply(spec.table, 1, form)


# convert to Mg./hectare
spec.table$biom <- biomass * spec.table$density / 1000

# do some removing of places with only one tree
one_tree <- c(which(!species[,1] %in% "No tree" & species[,2] %in% "No tree"),
              which(species[,1] %in% "No tree" & !species[,2] %in% "No tree"))

spec.table <- subset(spec.table, !point %in% one_tree)
                          
#spec.table <- spec.table[,2:14]
colnames(spec.table)[1:2] <- c("x", "y")# rename grid cell x and y colnames
write.csv(spec.table, 
        file = paste0('outputs/density_biomass_pointwise.ests_inilmi','_v', 
                      version, 
                      '.csv'), row.names = FALSE)

# there are some excessivly high estimates of density + biomss, so we convert these to the 99th percentile value:
                         
pre.quantile <- spec.table

#take the 99 percentile of these, since density blows up in some places
nine.nine.pct <- apply(spec.table[,c("density", "basal", "diams", "dists", "biom")], 2, quantile, probs = 0.995, na.rm=TRUE)

 
nine.five.pct <- apply(spec.table[,c("density", "basal", "diams", "dists", "biom")], 2, quantile, probs = 0.95, na.rm=TRUE)



# assign all points greater than the 99th percentile to 99th percentile values
spec.table$density[spec.table$density > nine.nine.pct['density']] <- nine.nine.pct['density']
spec.table$basal[spec.table$basal > nine.nine.pct['basal']] <- nine.nine.pct['basal']
spec.table  <- spec.table[!is.na(spec.table$density), ]

write.csv(spec.table, file=paste0('outputs/biomass_no_na_pointwise.ests_inilmi','_v',version, '.csv'), row.names = FALSE)

}
