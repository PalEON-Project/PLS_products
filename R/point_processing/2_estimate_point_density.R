## to be based on KAH 01b_calculate_full_density.r (and 02_calculate_density.R?)
## calculate point-level density estimates

library(plyr)
library(reshape2)
library(raster)
library(tidyr)
version <- "1.7-5"

library(plyr)
library(reshape2)
library(raster)
library(tidyr)
version <- "1.7-5"

#-----------------------load data------------------------------------------------

#read in final.data from the step_one_clean_IN.r script:
final.data <- read.csv(paste0("outputs/ndilin_pls_for_density_v",version,".csv"), stringsAsFactors = FALSE)
final.data <- final.data[,1:24]

# read in final data from michigan
final.data.mi <- read.csv(paste0("data/lower_mi_final_data.csv"), stringsAsFactors = FALSE)

# read in UMW data:
final.data.umw <- read.csv("data/outputs/Point_Data_From_Goringetal16_used_data_alb.csv", stringsAsFactors = FALSE)

# Combine all the data:

final.data <- rbind(final.data, final.data.mi, final.data.umw)


# read in charlies correction factor matching table:

corr.vals <- read.csv('data/charlie_corrections_full_midwest_mi_all.csv') # this file needs to be versioned

# match the correction factors to the combined dataset
match.vec <- apply(corr.vals[,c("state", "year", "corner", "sectioncorner")], 1, paste, collapse = '')
to.match <- apply(data.frame(final.data$state, final.data$surveyyear, final.data$corner, final.data$sectioncorner, stringsAsFactors = FALSE), 1, paste, collapse = '')

correction.factor <- corr.vals[match(to.match, match.vec),]


# also get species dataframe!
species <- final.data[,13:16]


#------------------------Estimate Tree Density-----------------------------------
## Morisita estimates for  densities and basal area with charlies correction factors

# morisita density estimator from Simon Goring's Witness tree code here: WitnessTrees/R/paper/R/misc.functionsv1.4.R

morisita <- function(processed.data, correction.factor = NULL, veil=FALSE) {
      #  Function to calculate stem density using the morista function.  The input 
      #  is 'processed.data', which should be the file 'used.data'.  'correction 
      #  factor' is the modified Cottam Correction factor determined in 
      #  'load.estimate.correction.R' using a generalized linear model (with a Gamma
      #  distribution).


      azim <- processed.data[,c('az1', 'az2', 'az3', 'az4')]
      diam <- processed.data[,c('diam1', 'diam2', 'diam3', 'diam4')]
      dist <- processed.data[,c('dist1', 'dist2', 'dist3', 'dist4')]
      spec <- processed.data[,c('species1', 'species2', 'species3', 'species4')]
      corn <- processed.data[,'corner']


      if(veil){
        diam[diam < 8] <- NA #double check that this is corrected veil line, and do we want to use it?
      }

      diam[diam == 0 & !spec == 'No tree'] <- NA

      #m.diam <- diam/100 #diameters are in cm already
      m.diam <- (diam*2.54 ) / 100 # convert diameter from in to meters

      dist <- floor(apply(dist, 2, function(x)as.numeric(as.character(x))))
      azim <- floor(apply(azim, 2, function(x)as.numeric(as.character(x))))

      #  This tells us how many quadrats are used.  I'd prefer to use all points
      #  where samples are drawn from two quadrats, but in some cases it seems that
      #  there are NAs in the data.
      #  If a point has recorded azimuths we state that they must be in two different
      #  quadrats:

      two.quads <- apply(azim[,1:2], 1, function(x) sum(!is.na(unique(floor(x/90)))))

      #  There are 10,155 points for which the first two trees were sampled in the
      #  same quadrat.  In general these are randomly distributed, but interestingly
      #  there's a big clump of them in Wisconsin.  Even so, there are lots of other
      #  points around.  We can accept that these points are categorically wrong.
      #  
      #  
      #  sum((two.quads == 1 & !(is.na(azim[,1]) | is.na(azim[,2]))))


      #  we need to change:

      two.quads[((two.quads < 2 & (is.na(azim[,1]) | is.na(azim[,2]))) &
                   !(is.na(dist[,1]) | is.na(dist[,2])))] <- 2

      #  Exclusions include:
      #  Plots with a tree as plot center:
      two.quads[dist[,1] == 0] <- 0

      #  Plots where one of the trees has no measured diameter:
      two.quads[is.na(diam[,1]) | is.na(diam[,2])] <- 0

      #  Plots where a distance to tree is missing:
      two.quads[is.na(dist[,1]) | is.na(dist[,2])] <- 0
      #  This is the same as k in Charlie's spreadsheet:

      q <- two.quads

      #  Tree dist is measured in links in the dataset, I am converting to
      #  meters and adding one half a dimater (in cm), on Charlie's advice.



         m.dist <- dist * 0.201168 + 0.5 * m.diam # convert distances from chains (links) to meters
         #m.dist <- dist + 0.5 * m.diam 
      #  rsum is the sum of the squared radii, in cases where there are two trees in
      #  the same quadrant I'm going to drop the site, as I will with any corner 
      #  with only one tree since the morista density estimator can't calculate
      #  density with less than two trees, and requires two quadrats.

      #  I'm going to let the NAs stand in this instance.
      rsum <- rowSums((m.dist[,1:2])^2, na.rm=T)
      #rmax <- max(m.dist[,1:2], na.rm=T)
      #rmax<- apply(m.dist[,1:2], 1, max, na.rm= T)

      #  A set of conditions to be met for the rsum to be valid:
      rsum[rowSums(is.na(m.dist[,1:2])) == 2 |  q < 2 | rsum == 0 | rowSums(m.dist[,1:2], na.rm=T) < 0.6035] <- NA
      #rmax[rowSums(is.na(m.dist[,1:2])) == 2 |  q < 2 | rmax == 0 | rowSums(m.dist[,1:2], na.rm=T) < 0.6035] <- NA

      #  From the formula,
      #  lambda = kappa * theta * (q - 1)/(pi * n) * (q / sum_(1:q)(r^2))
      #  here, n is equal to 1.
      #  units are in stems / m^2

      #charilies morisita.est have separate corrections for each type of corner
      morisita.est <- ((q - 1) / (pi * 1)) * (2 / rsum) *
        correction.factor$kappa  * correction.factor$theta* correction.factor$zeta * correction.factor$phi

      morisita.est[q < 2] <- NA

      #  Now they're in stems / hectare
      morisita.est <- morisita.est * 10000

      #  Basal area is the average diameter times the stem density.
      #  The stem density is measured in trees / ha.
      #met.rad <- (diam / 2)*2.54 / 100
      met.rad <- (diam / 2)/100
      basal.area <- morisita.est * rowSums(pi * met.rad^2, na.rm=TRUE)

      basal.area[ q < 2 ] <- NA



     # radius <- max(m.dist, na.rm = TRUE)
      #plot.area <- pi*rmax/2 # calculate average plot area
      return(list(morisita.est, basal.area))

}

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
                         township = final.data$Township)#,
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

