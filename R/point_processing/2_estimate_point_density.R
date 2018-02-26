## to be based on KAH 01b_calculate_full_density.r (and 02_calculate_density.R?)
## calculate point-level density estimates

library(plyr)
library(reshape2)
library(raster)
library(tidyr)
version <- "1.7-5"

#-----------------------load data------------------------------------------------

# read in pre-cleaned PLS point data from Indiana and Illinois
final.data <- read.csv(paste0("outputs/ndilin_pls_for_density_v",version,".csv"), stringsAsFactors = FALSE)

# read in the corrections for stem density for each grid cell:
correction.factor <- read.csv("data//correction_factors.csv", header = TRUE)
final.data <- final.data[,1:23] 

# read in final pre-cleaned data from southern michigan
final.data.mi <- read.csv(paste0("data/lower_mi_final_data.csv"), stringsAsFactors = FALSE)
final.data.mi <- final.data.mi[!names(final.data.mi) %in% c("cornertype", "NA.")] 

# TAKEN FROM SIMON's OUTPUT: read in final pre-cleaned data from the upper midwest:


# cogbill corrections for southern Michigan:
correction.factor.mi <- read.csv("data//MI_correction_factors.csv", header = TRUE)
correction.factor.mi <- correction.factor.mi[!names(correction.factor.mi) %in% c("X.1", "State", "Internal", "Section")]
colnames(correction.factor.mi) <- c("X","Pair", "kappa", "theta", "zeta", "phi")


# combine the southrn MI data and the INIL data: 
final.data <- rbind(final.data, final.data.mi)
correction.factor <- rbind(correction.factor, correction.factor.mi)# combine the southrn MI correction factors and the INIL correction factors: 


#------------------------Estimate Tree Density-----------------------------------
## Morisita estimates for indiana densities and basal area with charlies correction factors
# & no diameter veil
source('R/morisita.r') # morisita density estimator from Simon Goring's Witness Trees code

# morisita function calculates basal area and stem density
estimates <- morisita(final.data, correction.factor, veil = TRUE)

stem.density <- estimates[[1]]
basal.area <- estimates[[2]]

# there are some very high estimates of stem density 
summary(stem.density)
summary(basal.area)
zero.trees <- is.na(stem.density) 

#set stem.density = 0 where there are zero trees due to 'No tree' or 'Wet' or 'Water' as the species
stem.density[zero.trees] <- 0
basal.area[zero.trees] <- 0


# make into a data frame and export as csv
stem.density <- data.frame(stem.density, basal.area, final.data)
write.csv(stem.density, paste0('outputs/IN_IL_MI_densestimates_v',version,'.csv'))

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



