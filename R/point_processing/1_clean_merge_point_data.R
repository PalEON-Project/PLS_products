## to be based on KAH 01a_clean_merge_IN_IL.r 
## first step is to clean all the data for Southern MI, Uppermidwest, and Indiana + Illinois separately, get correction factors for all the data, then join together
## will work on combining the all the data before estimating correction factors, to make the correction factor generation more intuitive
## merge point data from UMW, IL, IN, SO MI

library(sp)
library(spdep)
library(rgdal)
library(raster)
library(ggplot2)
library(Rcpp)


# ----------------------------------DATA CLEANING: IN + IL --------------------------------------------------
#---------------------read in data and clean up the column names------------------------
version <- "1.7-5" # version using 1.8 IL data and 1.7 IN data

# Read in the data
ind <- read.csv("data/ndinpls_v1.7.csv", stringsAsFactors = FALSE) # version 1.6-1 


# Read in the il data on version v1.8
il <- read.csv("data/ndilpls_v1.8.csv", stringsAsFactors = FALSE) # version 1.6



#-------------------------Data cleaning--------------------------------------------------

# we can't use datapoints listed as no data, or those that are missing data in key variables
# converting all No data trees to NA's :
ind[ind$L3_tree1 %in% 'No data',] <- NA
ind[ind$L3_tree2 %in% 'No data',] <- NA

il[il$L3_tree1 %in% 'No data',] <- NA
il[il$L3_tree2 %in% 'No data',] <- NA


# remove all instances of no data for the l3tree 1
ind <- ind[!is.na(ind$L3_tree1),]
il <- il[!is.na(il$L3_tree1),]


#  IN distances are in chains
ind$DIST1 <- as.numeric(ind$chainstree)
ind$DIST2 <- as.numeric(ind$chainstree2)
ind$DIST3 <- as.numeric(ind$chainstree3)
ind$DIST4 <- as.numeric(ind$chainstree4)

# Il distances in chains to tree
il$DIST1 <- as.numeric(il$chainstree)
il$DIST2 <- as.numeric(il$chainstree2)
il$DIST3 <- as.numeric(il$chainstree3)
il$DIST4 <- as.numeric(il$chainstree4)

# make sure township names have the state in front of them:
ind$twp <- c(paste('IN', as.character(ind$TRP)))

il$twp_il <- c(paste('IL', as.character(il$TRP)))


## for the getAngle function to work later, we need a 4 character Azimuth from the tree bearings and bearing direction
ind$bearings1 <- c(paste0(as.character(ind$bearing),  as.character(ind$bearingdir)))
ind$bearings2 <- c(paste0(as.character(ind$bearing2),  as.character(ind$bearingdir2)))
ind$bearings3 <- c(paste0(as.character(ind$bearing3),  as.character(ind$bearingdir3)))
ind$bearings4 <- c(paste0(as.character(ind$bearing4),  as.character(ind$bearingdir4)))

il$bearings1 <- c(paste0(as.character(il$bearing),  as.character(il$bearingdir)))
il$bearings2 <- c(paste0(as.character(il$bearing2),  as.character(il$bearingdir2)))
il$bearings3 <- c(paste0(as.character(il$bearing3),  as.character(il$bearingdir3)))
il$bearings4 <- c(paste0(as.character(il$bearing4),  as.character(il$bearingdir4)))

il$state <-'IL'
ind$state <-'IN'

#create and rename columns to match up columns from indiana and illinois
il$twp <- il$TRP

il$DIST4 <- NA

keeps <- c("x","y","twp","year","L3_tree1", "L3_tree2", "L3_tree3", "L3_tree4", "bearing", 
  "bearing2", "bearing3", "bearing4","degrees", "degrees2", "degrees3","degrees4", "DIST1", "DIST2", "DIST3", "DIST4",
  "diameter", "diameter2", "diameter3", "diameter4", "cornerid", "typecorner","state")

ind.data <- ind[keeps] 
il.data <- il[keeps]

#  The merged dataset is called inil
inil <- rbind(data.frame(ind.data), data.frame(il.data))
inil <-data.frame(inil, stringsAsFactors = FALSE)


#  There are a set of 99999 values for distances which I assume are meant to be NAs. 

inil[inil == 88888 ] <- NA
inil[inil == 99999 ] <- NA
  
inil$bearing[inil$bearing == ''] <- NA     
inil$bearing2[inil$bearing2 == ''] <- NA
inil$bearing3[inil$bearing3 == ''] <- NA     
inil$bearing4[inil$bearing4 == ''] <- NA
inil$year[inil$year == 99999] <- NA # our correction factors are by year, so we need the year


# There are some points in Illinois where distances are listed as 0, but they are "Water" or "wet" or "No tree"
# Here we change these distnces to 'NA'

summary(inil[inil$L3_tree1 %in% c('No tree', 'Water', 'Wet') | inil$L3_tree2 %in% c('No tree', 'Water', 'Wet'),])
zero.trees <-(inil$L3_tree1 %in% c('No tree', 'Water', 'Wet') | inil$L3_tree2 %in% c('No tree', 'Water', 'Wet'))

inil[zero.trees, c("DIST1", "DIST2", "DIST3")] <- NA
inil[zero.trees, c('diameter', 'diameter2', "diameter3")] <- NA

# now kill missing cells:
inil <- inil[!is.na(inil$y),]
inil <- inil[!is.na(inil$x),]

# create a survey year variable that coresponds to survey year correction factors
year <- ifelse(inil$year >= 1825, '1825+',
                    ifelse(inil$year < 1825, '< 1825',"ALL"))

inil$surveyyear <- year

inil <- data.frame(inil)

# ----------------------------reorganizing INIL data -------------------------------------

# create data frames for diameters, distances, bearings and degrees
#diameters convertedto cm
diams <-  cbind(as.numeric(inil$diameter), 
                as.numeric(inil$diameter2), 
                as.numeric(inil$diameter3), 
                as.numeric(inil$diameter4))

#distances converted to meters
dists <-  cbind(as.numeric(inil$DIST1), 
                as.numeric(inil$DIST2), 
                as.numeric(inil$DIST3), 
                as.numeric(inil$DIST4))

bearings <- cbind(as.character(inil$bearing), 
                  as.character(inil$bearing2),
                  as.character(inil$bearing3),
                  as.character(inil$bearing4))

degrees <- cbind(as.numeric(inil$degrees), 
                 as.numeric(inil$degrees2),
                 as.numeric(inil$degrees3),
                 as.numeric(inil$degrees4))



#--------------------geting azimuths from distance and direction-----------------

#  Use Simon's getAngle function to find the azimuth 
#  getAngle converts the four character azimuth (e.g. N43E) to a numeric, 360
#  degree angle.  It also has to deal with a number of special cases.
#  The code for getAngles is a bit scuzzy, but it leaves only 231 azimuths 
#  untranslated, this is a manageable number.
get_angle_IN <- function(bearings, degrees, dists) {
  #  This function is used in 'step.one.clean.bind_v1.1.R', it converts the
  #  text azimuth strings to numeric, 360 degree values.
  #  This is the vector that will store the values.
  angl <- degrees
  
  #  This is a special case, generally where the tree is plot center.
  angl[degrees == '0' & dists =='0'] <- 0
  
  #  This is a short function that takes cares of NAs in boolean functions, it's
  #  just a simple wrapper for the boolean function that sets the NA values in
  #  the vector to FALSE.
  fx.na <- function(x) { x[ is.na( x ) ] <- FALSE; x }
  
  #  Given the text azimuths in the dataset, return the quadrant values.
  #  This gives a boolean index of the quadrant
  
  north <- fx.na(bearings == 'NNA'| bearings == 'NAN' | bearings =='N'|bearings =='N99999'|bearings =='N88888')
  east <- fx.na(bearings == 'NAE' | bearings =="ENA" | bearings =='E'|bearings =='99999E'|bearings =='88888E')
  south <- fx.na(bearings == 'SNA' | bearings =="NAS"| bearings == 'S'|bearings =='S99999'|bearings =='S88888')
  west <- fx.na(bearings == 'NAW' | bearings =="WNA" | bearings == 'W'|bearings =='99999W'|bearings =='88888W')
  #north <- fx.na( regexpr('N', bearings) > 0 )
  #east  <- fx.na( regexpr('E', bearings) > 0 | bearings == 'EAST')
  #south <- fx.na( regexpr('S', bearings) > 0 | bearings == 'SOUTH')
  #west <-  fx.na( regexpr('W', bearings) > 0 | bearings == 'WEST')
  
  ne <- fx.na( (north & east) | bearings == 'NE')
  se <- fx.na( (south & east) | bearings == 'SE')
  sw <- fx.na( (south & west) | bearings == 'SW')
  nw <- fx.na( (north & west) | bearings == 'NW') 
  
  #  The cell is in a quadrant, regardless of which.
  quad <- ne | se | sw | nw
  
  #  Special case of the trees with a unidirectional direction.
  uni  <- (!quad) & !(north | south | east | west) 
  
  angl[ uni & north ] <- 0
  angl[ uni & south ] <- 180
  angl[ uni & east  ] <- 90 
  angl[ uni & west  ] <- 270
  
 
  ##########
  
  #  Another set of special cases:
 
  degrees <- apply(degrees, 2, as.numeric)
  angl[ north ] <- degrees[ north ]
  angl[ east ] <- 180 - degrees[ east ]
  angl[ south ] <- 180 + degrees[ south ]
  angl[ west ] <- 360 - degrees [ west ]
  
  return(angl)
  
}



azimuths <- get_angle_IN(bearings, degrees, dists)

#####  Cleaning Trees:  This is already done in the CSV file, but we should check this with jody's update teo the paleon conversion file
#      Changing tree codes to lumped names:
#spec.codes <- read.csv('data/input/relation_tables/fullpaleon_conversion_v0.3-3.csv', stringsAsFactor = FALSE)
#spec.codes <- subset(spec.codes, Domain %in% 'Upper Midwest')

#lumped <- data.frame(abbr = as.character(spec.codes$Level.1),
 #                    lump = as.character(spec.codes$Level.3a))

species.old <- data.frame(as.character(inil$L3_tree1), 
                          as.character(inil$L3_tree2), 
                          as.character(inil$L3_tree3), 
                          as.character(inil$L3_tree4), stringsAsFactors = FALSE)
species <- species.old #since species is already converted
#species <- t(apply(species.old, 1, 
 #                  function(x) lumped[match(tolower(x), tolower(lumped[,1])), 2]))


#  Now we assign species that don't fit to the 'No tree' category.
species[is.na(species)] <- 'No tree'

#  Here there needs to be a check, comparing species.old against species.
test.table <- table(unlist(species.old), unlist(species), useNA='always')
write.csv(test.table, 'data/outputs/clean.bind.test.csv')

######
#  Some annoying things that need to be done:
#  First, there are some points where the first tree has a distance of zero
#  since it is the plot center.  
#  In these cases, the first azimuth is sometimes listed in a strange way, either
#  it's an NA (obviously) or it's a strange value.  In any case, we need to
#  ensure the azimuth is something recognized, I set it to 0.  It doesn't really
#  matter though.

treed.center <- (dists[,1] == 0 & !is.na(azimuths[,1]) & diams[,1] > 0)
treed.center[is.na(treed.center)] <- FALSE

azimuths[treed.center,1] <- 0 #assign azimuth to 0

#  Another special case, two trees of distance 1. 
dists[rowSums(dists == 1, na.rm=T) > 1, ] <- rep(NA, 4)



#--------------Reorder the tree number by distance to the point-----------------

#  At this point we need to make sure that the species are ordered by distance
#  so that trees one and two are actually the closest two trees.


sp.levels <- levels(factor(unlist(species)))

species.num <- t(apply(species, 1, function(x) match(x, sp.levels)))

usable.data <- data.frame(diams,
                          dists,
                          species.num,
                          azimuths,
                          stringsAsFactors = FALSE)

rank.fun <- function(x){
  #  This is the function we use to re-order the points so that the closest is first, based on distances.
  
  test.dists <- as.vector(x[c(5:8)])
  
  ranker <- order(test.dists, na.last=TRUE) #order function sorts the distances
  
  return(x[c(ranker, ranker+4, ranker + 8, ranker + 12)])
}

colnames(usable.data) <- c(paste('diam', 1:4, sep =''),
                           paste('dist', 1:4, sep = ''), 
                           paste('species', 1:4, sep = ''),
                           paste('az', 1:4, sep = ''))

ranked.data <- matrix(nrow = nrow(usable.data),
                      ncol = ncol(usable.data))

for(i in 1:nrow(ranked.data)){
  if( sum(is.na(ranked.data[i,5:8]))<2 ){
    ranked.data[i,] <- unlist(usable.data[i,]) # if there is only 1 tree, just use the usable data as is
  } else{
    ranked.data[i,] <- unlist(rank.fun(usable.data[i,])) #if there is more than 1 tree, rank.fun sorts by distance
  }
  if(i%%6500 == 0)cat(':)') #prints a ':)' for each 6500 rows
}

ranked.data <- t(apply(usable.data, 1, rank.fun)) # need to drop 'id'
ranked.data <-data.frame(ranked.data)

# Convert species from numeric codes back into text
species <- data.frame(species1 = sp.levels[ranked.data[, 9]],
                      species2 = sp.levels[ranked.data[,10]],
                      species3 = sp.levels[ranked.data[,11]],
                      species4 = sp.levels[ranked.data[,12]],
                      stringsAsFactors=FALSE)
                       
 
#----------Getting correction factors----------------------



#  Indiana data has same correction factors for the whole state
# correction factors vary depending on which type of corner you are at
                       
extsec <- c(100100,200100, 300100, 400100, 500100, 600100, 700100,
            100200, 100300,100400,100500, 100600, 100700, 
            200700, 300700, 400700, 500700, 600700, 700700, 
            700100, 700200, 700300, 700400, 700500, 700600)
extqtr <- c(100140, 100240, 100340, 100440, 100540, 100640,
            140100, 240100, 340100, 440100, 540100, 640100, 
            140700, 240700, 340700, 440700, 540700, 640700, 
            700140, 700240, 700340, 700440, 700540, 700640)
intsec <- c(200200, 300200, 400200, 500200, 600200,
            200300, 300300, 400300, 500300, 600300,
            200400, 300400, 400400, 500400, 600400,
            200500, 300500, 400500, 500500, 600500,
            200600, 300600, 400600, 500600, 600600)
intqtr <- c(140200, 240200, 340200, 440200, 540200, 640200,
            140300, 240300, 340300, 440300, 540300, 640300,
            140400, 240400, 340400, 440400, 540400, 640400,
            140500, 240500, 340500, 440500, 540500, 640500,
            140600, 240600, 340600, 440600, 540600, 640600,
            200140, 300140, 400140, 500140, 600140,
            200240, 300240, 400240, 500240, 600240,
            200340, 300340, 400340, 500340, 600340,
            200440, 300440, 400440, 500440, 600440,
            200540, 300540, 400540, 500540, 600540,
            200640, 300640, 400640, 500640, 600640)
corner <- rep('NA', length(inil$cornerid))

corner <- ifelse(inil$cornerid %in% intsec, 'intsec',
       ifelse(inil$cornerid %in% intqtr, 'intqtr',
              ifelse(inil$cornerid %in% extsec, 'extsec',
                     ifelse(inil$cornerid %in% extqtr,  'extqtr',
                            ifelse(inil$typecorner == "(1/4) Section", "intqtr",
                                   ifelse(inil$typecorner == "Section", "intsec", "extsec"))))))
       

inil$cornertype <- paste0(corner, inil$state)


#These are the columns for the final dataset.
internal <- ifelse(!inil$cornerid %in% c("extsec", "extqtr", "external"), 'internal', 'external')
#trees    <- ifelse(plot.trees == 2, 'P', '2NQ')
section  <- ifelse(inil$typecorner %in% "Section", 'section', 'quarter-section')
point <- rep("P", length(inil$typecorner))

final.data <- data.frame(inil$x,
                         inil$y,
                         inil$twp,
                         as.character(inil$state),
                         ranked.data[,1:8],
                         species[,1:4],
                         ranked.data[,13:16], 
                         internal,
                         section,
                         inil$surveyyear,
                         point,
                         stringsAsFactors = FALSE)

colnames(final.data) <- c('PointX','PointY', 'Township','state',
                          paste('diam',    1:4, sep =''),
                          paste('dist',    1:4, sep = ''), 
                          paste('species', 1:4, sep = ''),
                          paste('az',      1:4, sep = ''), 'corner', "sectioncorner",'surveyyear', "point")

                          
write.csv(full.final, paste0("outputs/ndilin_pls_for_density_v",version,".csv"))


# ----------------------------------DATA CLEANING: SOUTHERN MI --------------------------------------------------

mich <- read.csv("data/southernmi_projected_v1/southernMI_projected_v1.0.csv", stringsAsFactors = FALSE)


not.no.tree <- !(!is.na(mich$L3_tree1) & is.na(mich$species1))
no.tree     <- is.na(mich$species1)
#mich <- mich[not.no.tree & !no.tree,]

#  Character vectors are read into R as factors, to merge them they need to
#  be converted to character strings first and then bound together.  To ensure
#  township and ranges are unique we add a state abbreviation to the front of
#  the twonship name.
twp <- c(#paste('mn', as.character(minn$TWP)), 
         #paste('wi', as.character(wisc$TOWNSHIP)), 
         paste('mi', as.character(mich$town)))
rng <- c(#as.character(minn$RNG), 
         #paste(as.character(wisc$RANGE),as.character(wisc$RANGDIR), sep=''), 
         as.character(mich$range))

#  The merged dataset is called nwmw, Minnesota comes first, then Wisconsin.
#nwmw <- rbind(minn[,c(8, 10:25)], wisc[,c(5, 13:28)], mich[,c(36, 13:28)])
nwmw <- mich
nwmw$twp <- twp
nwmw$rng <- rng

 #  There are a set of 9999 values for distances which I assume are meant to be NAs.  Also, there are a set of points where
#  the distance to the tree is 1 or 2 feet.  They cause really big density estimates!
nwmw [ nwmw == '9999'] <- NA
nwmw [ nwmw == '8888'] <- NA
nwmw [ nwmw == '_'] <- NA       # Except those that have already been assigned to 'QQ'
nwmw [ nwmw == '99999'] <- NA
nwmw [ nwmw == '999999'] <- NA
nwmw [ nwmw == '6666'] <- NA
nwmw [ nwmw == '999'] <- NA

nwmw[(is.na(nwmw$species1) & nwmw$diam1 > 0) | (is.na(nwmw$species2) & nwmw$diam2>0),] <- rep(NA, ncol(nwmw))  #  removes four records with no identified trees, but identified diameters

diams <-  cbind(as.numeric(nwmw$diam1), 
                as.numeric(nwmw$diam2), 
                as.numeric(nwmw$diam3), 
                as.numeric(nwmw$diam4))

dists <-  cbind(as.numeric(nwmw$dist1), 
                as.numeric(nwmw$dist2), 
                as.numeric(nwmw$dist3), 
                as.numeric(nwmw$dist4))

# mi azimuths are from 0 to 60
azimuths <- cbind(as.numeric(nwmw$az1_360), 
                  as.numeric(nwmw$az2_360),
                  as.numeric(nwmw$az3_360),
                  as.numeric(nwmw$az4_360))

colnames(azimuths) <- c("az1", "az2","az3", "az3")


species <- cbind(as.character(nwmw$L3_tree1), 
                 as.character(nwmw$L3_tree2), 
                 as.character(nwmw$L3_tree3), 
                 as.character(nwmw$L3_tree4))


species[species %in% ''] <- 'No tree'

#  Now we assign species that don't fit to the 'No tree' category.
species[is.na(species)] <- 'No tree'


######
#  Some annoying things that need to be done:
#  First, there are some points where the first tree has a distance of zero
#  since it is the plot center.  
#  In these cases, the first azimuth is sometimes listed in a strange way, either
#  it's an NA (obviously) or it's a strange value.  In any case, we need to
#  ensure the azimuth is something recognized, I set it to 0.  It doesn't really
#  matter though.

treed.center <- (dists[,1] == 0 & !is.na(azimuths[,1]) & diams[,1] > 0)
treed.center[is.na(treed.center)] <- FALSE

azimuths[treed.center,1] <- 0

#  Another special case, two trees of distance 1.  What's up with that?!
dists[rowSums(dists == 1, na.rm=T) > 1, ] <- NA#,rep(NA, 4)

#  When the object is NA, or the species is not a tree (NonTree or Water), set
#  the distance to NA.
dists[is.na(species) | species %in% c('No tree', 'Water', 'Missing')] <- NA


#  At this point we need to make sure that the species are ordered by distance
#  so that trees one and two are actually the closest two trees.


sp.levels <- levels(factor(unlist(species)))

species.num <- t(apply(species, 1, function(x) match(x, sp.levels)))

usable.data <- data.frame(diams,
                          dists,
                          species.num,
                          azimuths,
                          stringsAsFactors = FALSE)

rank.fun <- function(x){
  #  This is the function we use to re-order the points so that the closest is first, based on distances.
  
  test.dists <- as.vector(x[c(5:8)])
  
  ranker <- order(test.dists, na.last=TRUE)
  
  return(x[c(ranker, ranker+4, ranker + 8, ranker + 12)])
}

colnames(usable.data) <- c(paste('diam', 1:4, sep =''),
                         paste('dist', 1:4, sep = ''), 
                         paste('species', 1:4, sep = ''),
                         paste('az', 1:4, sep = ''))

ranked.data <- matrix(nrow = nrow(usable.data),
                      ncol = ncol(usable.data))

for(i in 1:nrow(ranked.data)){
  if( sum(is.na(ranked.data[i,5:8]))<2 ){
    ranked.data[i,] <- unlist(usable.data[i,])
  } else{
    ranked.data[i,] <- unlist(rank.fun(usable.data[i,]))
  }
  if(i%%6500 == 0)cat('.')
}

ranked.data <- t(apply(usable.data, 1, rank.fun)) # need to drop 'id'
colnames(ranked.data) <- c(paste('diam', 1:4, sep =''),
                           paste('dist', 1:4, sep = ''), 
                           paste('species', 1:4, sep = ''),
                           paste('az', 1:4, sep = ''))

species <- data.frame(species1 = sp.levels[as.numeric(ranked.data[, 9])],
                      species2 = sp.levels[as.numeric(ranked.data[,10])],
                      species3 = sp.levels[as.numeric(ranked.data[,11])],
                      species4 = sp.levels[as.numeric(ranked.data[,12])])

year <- rep("NA", length(species$species1))
state <- data.frame(state = rep("MI", length(species$species1)))

corner <- mich$sec_corner

state <- rep("Michigan", length(species$species1))
year <- rep("NA", length(species$species1))

# assign year as the SW or SE of michigan (this is just for correction factors)
year[state == 'Michigan' & nwmw$twnrng %like% "W"] <- 'SW'
year[state == 'Michigan' & nwmw$twnrng %like% "E"] <- 'SE'

plot.trees <- rowSums(!(species == 'Water' | species == 'No tree'), na.rm = TRUE)

point.no <- as.character( mich$sec_corner)

#  So there are a set of classes here, we can match them all up:

internal <- ifelse(!point.no %in% "Extsec", 'internal', 'external')
trees    <- ifelse(plot.trees == 2, 'P', '2nQ')
section  <- ifelse(mich$sec_corner %in% "section", 'section', 'quarter-section')


#  These are the columns for the final dataset.

final.data <- data.frame(nwmw$point_x,
                         nwmw$point_y,
                        nwmw$twnrng,
                        state,
                        ranked.data[,1:8],
                        species,
                        ranked.data[,13:16],
                        internal,
                        section,
                        year ,
                        trees,
                        stringsAsFactors = FALSE)

colnames(final.data) <- c('PointX','PointY', 'Township',"state",
                          paste('diam',    1:4, sep =''),
                          paste('dist',    1:4, sep = ''), 
                          paste('species', 1:4, sep = ''),
                          paste('az',      1:4, sep = ''), 'corner',"sectioncorner",'surveyyear', 'point')


final.data$az1[final.data$az1 <= 0 ] <- NA
final.data$az2[final.data$az2 <= 0] <- NA
final.data$az3[final.data$az3 <= 0] <- NA
final.data$az4[final.data$az4 <= 0] <- NA

summary(final.data)
final.data <- final.data[!is.na(final.data$PointX),]
                       
#write to a csv:
write.csv(final.data, "data/lower_mi_final_data.csv")


                       
# ----------------------------------DATA CLEANING: UMW -------------------------------------------------------
 # data cleaning modified from Simon's witness tree code (https://github.com/PalEON-Project/WitnessTrees/blob/master/R/process_raw/step.one.clean.bind_v1.4.R)
 # correction factors for UMW are also generated in the witness tree code

#  Binding and cleaning the Wisconsin and Minnesota Public Lands Surveys, data is
#  sourced from the Mladenoff Lab at the University of Wisconsin.  The lab has
#  granted us permission to use the data, but not permission to distribute the
#  original datasets.  A version of the Minnesota data can be obtained from 
#  http://deli.dnr.state.mn.us/metadata/pveg_btreept3.html
#  The wisconsin data may be obtained by contacting David Mladenoff at:
#  mladenoff@wisc.edu
#
#  This file opens the Wisconsin and Minnesota datasets, renames the columns of
#  the Minnesota shapefile to match those of the Wisconsin dataset, and then
#  binds a number of columns from both datasets together (but not the complete
#  set of columns, since there are some columns unique to each dataset.  The
#  The ultimate dataset has the following columns:
#  Point:  PLS point number
#  Township:  Township line
#  Range:  Range line
#  diam (1 through 4):  Bearing tree diameter
#  dist  (1 through 4):  Distance to bearing tree
#  species_level1  (1 through 4):  Species as given in the Public Land Surveys for the bearing tree'
#  species_level3a  (1 through 4):  Taxa of the bearing tree based on the conversion file 
#                           'level0_to_level3a_v0.4-7.csv'
#  az (1 through 4): Azimuth to the bearing trees.

#  Data:  The following files are distributed with this code:
#  Maps/glo_corn_ex.shp : This dataset represents a subset of the full Wisconsin data,
#  and uses only 1% of the original dataset.
#  Maps/Minnesota_ex.shp : 5% of the original Minnesota dataset.

library(sp)
library(spdep)
library(rgdal)
library(raster)

wisc <- readOGR('data/wisc/glo_corn.shp','glo_corn')
minn <- readOGR('data/minn/Minnesota.shp', 'Minnesota')
mich <- readOGR('data/mich/michigan_filled/michigan_filled.shp', 'michigan_filled')

#  The files are in unique projections, this standardizes the projections to
#  a long-lat projection:

wisc <- spTransform(wisc, CRS('+proj=longlat +ellps=WGS84'))
minn <- spTransform(minn, CRS('+proj=longlat +ellps=WGS84'))
mich <- spTransform(mich, CRS('+proj=longlat +ellps=WGS84'))

#  The wisconsin Range is set as a single value, the 'E' and 'W' codes are in
#  RANGDIR.  Looking at the data it also looks like there are a few ranges
#  that are miscoded.
wisc$RANGDIR[wisc$RANGDIR %in% '2'] <- 'W'
wisc$RANGDIR[wisc$RANGDIR %in% '4'] <- 'E'

wisc$RANGDIR[wisc$RANGDIR %in% 'W' & coordinates(wisc)[,1] > 5e+05] <- 'E'

#  Michigan's point numbers are wrong in the dataset.  I'm not sure where the 
#  error arose from, but we need them to be able to assign section & quartersection
#  points.
mich$pnt <- as.numeric((substr(as.character(mich$RECNUM_C), 9,11)))

#  The index of column names in Minnesota becomes the same as the Wisconsin.
names(minn)[c(8, 10:25)] <- names(wisc)[c(5, 13:28)]
names(mich)[c(36, 13:28)] <- names(wisc)[c(5, 13:28)]

#  This step throws up four warnings, the warnings seem to come from 
#  the distance fields, which are stored as a factor for some reason.
minn$DIST1 <- as.numeric(levels(minn$DIST1)[minn$DIST1])
minn$DIST2 <- as.numeric(levels(minn$DIST2)[minn$DIST2])
minn$DIST3 <- as.numeric(levels(minn$DIST3)[minn$DIST3])
minn$DIST4 <- as.numeric(levels(minn$DIST4)[minn$DIST4])

#  We have made a choice to say that all taxa labelled 'Beech' in Minnesota are likely
#  Bluebeech, or, in our dataset, Ironwood.
minn@data[minn@data == 'BE'] <- 'IR'

#  We want the Minnesota data to reflect water in the same way that the Wisconsin data does.
#  Almendinger references the following land cover codes that are likely to have water:
#  'A' - Creek (unlikely to be exclusively water)
#  'M' - Marsh
#  'S' - Swamp
#  'L' - Lake

minn$SP1 <- as.character(minn$SP1)
minn$SP2 <- as.character(minn$SP2)
minn$SP1[minn$VEGTYPE %in% c('L', 'M', 'S', 'R', 'A') & minn$SP1 %in% '_'] <- 'QQ'
minn$SP2[minn$VEGTYPE %in% c('L', 'M', 'S', 'R', 'A') & minn$SP2 %in% '_'] <- 'QQ'

#  There are also some weird Michigan points:
#  1.  Michigan has a set of points with NA as SPP1 but identifiable trees listed as
#      'tree'.  1549 of these are quartersection points, 45 are section points.  This
#      is clearly an artifact of the sampling method.  We remove these points.
#  2.  There are also 2909 'no tree points in Michigan.  Most of these points are quarter
#      section points, and there is clear grographic bias.  We assume these points are
#      early survey points and remove them entirely.
not.no.tree <- !(!is.na(mich$TREE) & is.na(mich$SP1))
no.tree     <- is.na(mich$SP1)
mich <- mich[not.no.tree & !no.tree,]

#  Character vectors are read into R as factors, to merge them they need to
#  be converted to character strings first and then bound together.  To ensure
#  township and ranges are unique we add a state abbreviation to the front of
#  the twonship name.
twp <- c(paste('mn', as.character(minn$TWP)), 
         paste('wi', as.character(wisc$TOWNSHIP)), 
         paste('mi', as.character(mich$twp)))
rng <- c(as.character(minn$RNG), 
         paste(as.character(wisc$RANGE),as.character(wisc$RANGDIR), sep=''), 
         as.character(mich$rng))

#  The merged dataset is called nwmw, Minnesota comes first, then Wisconsin.
nwmw <- rbind(minn[,c(8, 10:25)], wisc[,c(5, 13:28)], mich[,c(36, 13:28)]) #six invalid factor level warnings come up

nwmw$twp <- twp
nwmw$rng <- rng

#  There are a set of 9999 values for distances which I assume are meant to be NAs.  Also, there are a set of points where
#  the distance to the tree is 1 or 2 feet.  They cause really big density estimates!
nwmw@data [ nwmw@data == '9999'] <- NA
nwmw@data [ nwmw@data == '8888'] <- NA
nwmw@data [ nwmw@data == '_'] <- NA       # Except those that have already been assigned to 'QQ'
nwmw@data [ nwmw@data == '99999'] <- NA
nwmw@data [ nwmw@data == '999999'] <- NA
nwmw@data [ nwmw@data == '6666'] <- NA
nwmw@data [ nwmw@data == '999'] <- NA

# There is some cleaning to do.  A bit frustrating.  We can't confirm the diameters of
#  a number of points, although we hope to at some point in the future:
#  No stem density removals, none of the plots look like they have 'weird' points.
#  Basal area removals:
nwmw@data[which(as.numeric(nwmw$DIAM1) >100),] <- rep(NA, ncol(nwmw))  #  removes 19 trees with reported diameters over 250cm.
nwmw@data[which(as.numeric(nwmw$DIAM2) >100),] <- rep(NA, ncol(nwmw))  #  removes an additional 14 trees.
nwmw@data[(is.na(nwmw$SP1) & nwmw$DIAM1>0) | (is.na(nwmw$SP2) & nwmw$DIAM2>0),] <- rep(NA, ncol(nwmw))  #  removes four records with no identified trees, but identified diameters

diams <-  cbind(as.numeric(nwmw$DIAM1), 
                as.numeric(nwmw$DIAM2), 
                as.numeric(nwmw$DIAM3), 
                as.numeric(nwmw$DIAM4))

dists <-  cbind(as.numeric(nwmw$DIST1), 
                as.numeric(nwmw$DIST2), 
                as.numeric(nwmw$DIST3), 
                as.numeric(nwmw$DIST4))

azimuths <- cbind(as.character(nwmw$AZ1), 
                  as.character(nwmw$AZ2),
                  as.character(nwmw$AZ3),
                  as.character(nwmw$AZ4))

#  getAngle converts the four character azimuth (e.g. N43E) to a numeric, 360
#  degree angle.  It also has to deal with a number of special cases.
#  The code for getAngles is a bit scuzzy, but it leaves only 231 azimuths 
#  untranslated, this is a manageable number.
source('R/point_processing/get_angle.R')
azimuths <- apply(azimuths, 2, get_angle)

#####  Cleaning Trees:  
#      Changing tree codes to lumped names:
spec.codes <- read.csv('data/level0_to_level3a_v0.4-7.csv', stringsAsFactor = FALSE)
spec.codes <- subset(spec.codes, domain %in% 'Upper Midwest')

lumped <- data.frame(abbr = as.character(spec.codes$level1),
                     lump = as.character(spec.codes$level3a))

species.old <- data.frame(as.character(nwmw$SP1), 
                          as.character(nwmw$SP2), 
                          as.character(nwmw$SP3), 
                          as.character(nwmw$SP4), stringsAsFactors = FALSE)

species <- t(apply(species.old, 1, 
                   function(x) lumped[match(tolower(x), tolower(lumped[,1])), 2]))

#  We need to indicate water and remove it.  There are 43495 cells with 'water'
#  indicated, and another 784 cells with 'missing' data.
#  when we limit these to the first two columns of the species table we get
#  a total of 25416 samples removed.

#  There are a set of dead taxa (DA, DB & cetera) that we exclude.  Only AM is
#  unknown at this point.  This excludes 213 trees.
species[species %in% ''] <- 'No tree'

#  Now we assign species that don't fit to the 'No tree' category.
species[is.na(species)] <- 'No tree'

#  Here Simon did a check comparing species.old against species.
#test.table <- table(unlist(species.old), unlist(species), useNA='always')
#write.csv(test.table, 'data/output/clean.bind.test.csv')

######
#  Some annoying things that need to be done:
#  First, there are some points where the first tree has a distance of zero
#  since it is the plot center.  
#  In these cases, the first azimuth is sometimes listed in a strange way, either
#  it's an NA (obviously) or it's a strange value.  In any case, we need to
#  ensure the azimuth is something recognized, I set it to 0.  It doesn't really
#  matter though.

treed.center <- (dists[,1] == 0 & !is.na(azimuths[,1]) & diams[,1] > 0)
treed.center[is.na(treed.center)] <- FALSE

azimuths[treed.center,1] <- 0

#  Another special case, two trees of distance 1.  What's up with that?!
dists[rowSums(dists == 1, na.rm=T) > 1, ] <- rep(NA, 4)

#  When the object is NA, or the species is not a tree (NonTree or Water), set
#  the distance to NA.
dists[is.na(species) | species %in% c('No tree', 'Water', 'Missing')] <- NA

#  Now we're looking at 36843 points with no usable data,
#  5022 points with only one tree
#  524 points with only two trees
#  59 points with three trees
#  380645 points with four trees

#  At this point we need to make sure that the species are ordered by distance
#  so that trees one and two are actually the closest two trees.

#  There's an annoying problem that has to do with having a character string in
#  the subsequent sort/order function, in that it converts everything to a
#  character string.  To fix it I change the string 'species' into a numeric
#  where each number is associated with a factor level.

sp.levels <- levels(factor(unlist(species)))

species.num <- t(apply(species, 1, function(x) match(x, sp.levels)))

usable.data <- data.frame(diams,
                          dists,
                          species.num,
                          azimuths,
                          stringsAsFactors = FALSE)


rank.fun <- function(x){
  #  This is the function we use to re-order the points so that the closest is first, based on distances.
  
  test.dists <- as.vector(x[c(5:8)])
  
  ranker <- order(test.dists, na.last=TRUE)
  
  return(x[c(ranker, ranker+4, ranker + 8, ranker + 12)])
}

colnames(usable.data) <- c(paste('diam', 1:4, sep =''),
                           paste('dist', 1:4, sep = ''), 
                           paste('species', 1:4, sep = ''),
                           paste('az', 1:4, sep = ''))
                           

ranked.data <- matrix(nrow = nrow(usable.data),
                      ncol = ncol(usable.data))

for(i in 1:nrow(ranked.data)){
  if( sum(is.na(ranked.data[i,5:8]))<2 ){
    ranked.data[i,] <- unlist(usable.data[i,])
  } else{
    ranked.data[i,] <- unlist(rank.fun(usable.data[i,]))
  }
  if(i%%6500 == 0)cat('.')
}

#ranked.data <- t(apply(usable.data, 1, rank.fun)) # need to drop 'id'

species <- data.frame(species1 = sp.levels[ranked.data[, 9]],
                      species2 = sp.levels[ranked.data[,10]],
                      species3 = sp.levels[ranked.data[,11]],
                      species4 = sp.levels[ranked.data[,12]])

#  We need to bin the year information so that we can use it to calculate
#  appropriate Cottam Correction factors.  The survey instructions for the PLS
#  change at a number of points during the sirveys in Wisconsin, but are
#  considered to be fixed by the time.
#  Some wisconsin samples don't have a year.  Look this up and figure out why.
#  It causes a problem with the Cottam correction factor.

mn_survey <- read.csv('data/minn/MN_Surveys.csv')
mn_survey$TOWN <- paste('T', formatC(mn_survey$TOWN, width=3, flag='0'), 'N', sep ='')
mn_survey$RANG <- paste('R', formatC(mn_survey$RANG, width=2, flag='0'), mn_survey$RDIR, sep ='')

get.minn.year <- function(x)which(paste(mn_survey$TOWN, mn_survey$RANG) == paste(x$TWP, x$RNG))

minn.year <- rep(NA, nrow(minn@data))
for(i in 1:nrow(minn@data)){
  if(is.na(minn.year[i])){
    minn.test <- get.minn.year(minn@data[i,])
    if(length(minn.test) == 1)minn.year[i] <- mn_survey$YEAR[minn.test]
  } 
}

wisc.year <- ifelse(wisc@data$YEAR_ > 1851, '1851+',
                    ifelse(wisc@data$YEAR_ > 1846, '1846-1851',
                           ifelse(wisc@data$YEAR_ > 1834, '1834-1846',
                                  ifelse(wisc@data$YEAR_ > 1832, '1832-1834','None'))))

#  Michigan has 47 instances where the year is '2', and 12 where the year is '9999'
#  The 9999s don't clean up because we're not using the nwmw data here:

mich.year <- as.numeric(as.character(mich@data$SURVYR))
mich.year[mich.year == '9999'] <- 2

mich.year <- ifelse(mich.year > 1851, '1851+',
                    ifelse(mich.year > 1846, '1846-1851',
                           ifelse(mich.year > 1834, '1834-1846',
                                  ifelse(mich.year > 1832, '1832-1834','None'))))
mich.year[is.na(mich.year)] <- 'None'

survey.year <- factor(c(minn.year, wisc.year, mich.year))

#  These are the columns for the final dataset.

final.data <- data.frame(nwmw$POINT,
                         nwmw$twp,
                         nwmw$rng,
                         ranked.data[,1:8],
                         species.old, 
                         species,
                         ranked.data[,13:16],
                         survey.year,
                         stringsAsFactors = FALSE)

colnames(final.data) <- c('Point', 'Township', 'Range',
                          paste('diam',    1:4, sep =''),
                          paste('dist',    1:4, sep = ''), 
                          paste('level1_',  1:4, sep = ''),
                          paste('level3a_', 1:4, sep = ''),
                          paste('az',      1:4, sep = ''), 'year')

#  Turn it into a SpatialPointsDataFrame and project into Great Lakes St.Lawrence Albers projection:
coordinates(final.data) <- coordinates(nwmw)
proj4string(final.data) <- proj4string(nwmw)
final.data <- spTransform(final.data, CRS('+proj=longlat +init=EPSG:3175'))

# now kill missing cells:
#final.data <- final.data[!final.data$species1 %in% c('Water', 'Missing'),] 
#final.data <- final.data[!final.data$species2 %in% c('Water', 'Missing'),]
#when Jody ran this all the data disappeared...



#  Write the data out as a shapefile.
writeOGR(final.data, 
         'data/output/uppermidwest_v1.shp', 
         'uppermidwestv1', 'ESRI Shapefile',
         overwrite_layer = TRUE, check_exists = TRUE)  
                       
                       
#write the data as a with the data and the coords
uppermidwest.coords = cbind(final.data@data, final.data@coords)
colnames(uppermidwest.coords) <- c('Point', 'Township', 'Range',
                          paste('diam',    1:4, sep =''),
                          paste('dist',    1:4, sep = ''), 
                          paste('level1_',  1:4, sep = ''),
                          paste('level3a_', 1:4, sep = ''),
                          paste('az',      1:4, sep = ''), 'year', 'x','y')

write.csv(uppermidwest.coords, 'data/output/uppermidwest.coords_v1.csv', row.names = FALSE)


                       
# kh: need to get the correction factors from UMW:
                       
                       
model.proj <- '+init=epsg:3175'
#  We use two different projection systems here.  This is the test to create the
#  base resolution.
if (model.proj == '+init=epsg:4326') {
  #  lat/long
  base.rast <- raster(xmn = -98.6, xmx = -66.1, ncols = 391,
                      ymn = 36.5,  ymx = 49.75, nrows = 160,
                      crs = '+init=epsg:4326')
  numbered.rast <- setValues(base.rast, 1:ncell(base.rast))
}

if (model.proj == '+init=epsg:3175') {
  base.rast <- raster(xmn = -71000, xmx = 2297000, ncols = 296,
                      ymn = 58000,  ymx = 1498000, nrows = 180,
                      crs = '+init=epsg:3175')
  numbered.rast <- setValues(base.rast, 1:ncell(base.rast))
}

if (model.proj == '+init=epsg:3175') {
  # These are narrower boundaries than the base raster for reasons associated with
  # the broader PalEON project.
  xylimits <- c(-100000, 1050000, 600000, 1600000)
}
if (model.proj == '+init=epsg:4326') {
  xylimits <- c(-98, -83, 42, 50)
}



used.data <- read.csv('data/output/uppermidwest.coords_v1.csv')

                       coordinates(used.data ) <- ~lon + lat
if (is.na(proj4string(used.data))) {
   #Just in case there happens to be no projection information associated with the data.
  proj4string(used.data) <- CRS('+proj=longlat +ellps=WGS84')
}

# converting to great lakes albers proj. and saving for jody
used.data.alb <- data.frame(spTransform(used.data, CRSobj = CRS("+init=epsg:3175")))



# need to rearrange data to match southern mi and in/il data:

colnames(used.data.alb)<- c("Point",      "Township" ,  "Range",      "diam1" ,     "diam2"  ,    "diam3" ,     "diam4",     
"dist1"     , "dist2",      "dist3"  ,    "dist4",      "species1",   "species2",   "species3",  
 "species4",   "az1"  ,      "az2"  ,      "az3"  ,      "az4",        "surveyyear",       "X",         
 "X.1"    ,    "optional",   "PointX"    ,    "PointY"   ,     "optional.1")

used.data.alb$state <- ifelse(substr(used.data@data$Township, 1, 2) == 'wi', 'Wisconsin',
                          ifelse(substr(used.data@data$Township, 1, 2) == 'mi', 'Michigan', 'Minnesota'))







#########################################################################
#  Clean the tree data:

diams <- used.data.alb[,c("diam1", "diam2", "diam3", "diam4")]
angles <- used.data.alb[,c("az1", "az2", "az3", "az4")]
dists <- floor(used.data.alb[,c("dist1", "dist2", "dist3", "dist4")])
species <- apply(used.data.alb[,c("species1", "species2", "species3", "species4")], 2, as.character)
species[is.na(species)] <- 'No tree'

#  Points within a township are either sections or quartersections.  This
#  is the list of points that are sections.  All others are quarter-sections.
sections <- c(2, 5, 8, 11, 14, 18, 21, 24, 27, 30,
              34, 37, 40, 43, 46, 50, 53, 56, 59, 62,
              66, 70, 74, 78, 82,
              87, 89, 91, 93, 95, 98, 100, 102, 104, 106, 108,
              109, 111, 113, 115, 117, 119, 122, 123, 124, 125, 126)

#  These are the points on the outside of each township.
external <- c(109:120, 97:108, 87, 89, 91, 93, 95, 122:126)

#  These correction values are derived empirically and are described in the supplementary material.
#  One issue right now is the lack of an empirical theta value.

#corr.vals <- read.csv('data/charlie_corrections_full_midwest_mi_all.csv')

correction <- data.frame(kappa = rep(NA, length(used.data.alb)),
                         theta = rep(NA, length(used.data.alb)),
                         zeta  = rep(NA, length(used.data.alb)),
                         phi   = rep(NA, length(used.data.alb)))

plot.trees <- rowSums(!(species == 'Water' | species == 'NonTree'), na.rm = TRUE)

point.no <- as.numeric(as.character(used.data.alb$Point))

#  So there are a set of classes here, we can match them all up:

internal <- ifelse(!point.no %in% external, 'internal', 'external')
trees    <- ifelse(plot.trees == 2, 'P', '2NQ')
cornersection  <- ifelse(point.no %in% sections, 'section', 'quartersection')
#corner <- ifelse(point.no %in% sections, 'section' & point.no %in% "external" , "Extsec","Intsec")
state    <- ifelse(substr(used.data.alb$Township, 1, 2) == 'wi', 'Wisconsin',
                   ifelse(substr(used.data.alb$Township, 1, 2) == 'mi', 'Michigan', 'Minnesota'))

used.data.alb$corner <- internal
used.data.alb$sectioncorner <- cornersection
used.data.alb$point <- trees

corr.year     <- as.character(used.data.alb$year)
corr.year[state == 'Michigan' ] <- "allN"
corr.year[state == 'Wisconsin' & corr.year %in% c('1832-1834', '1834-1846')] <- 1845
corr.year[state == 'Wisconsin' & !(corr.year %in% c('1832-1834', '1834-1846'))] <- 1907
corr.year[state == 'Minnesota' & (corr.year %in% (as.character(1847:1855)))] <- 1855
corr.year[state == 'Minnesota' & !(corr.year %in% (as.character(1847:1855)))] <- 1907
used.data.alb$surveyyear<- corr.year

used.data.alb <- used.data.alb[,c( "PointX", "PointY", "Township","state", "diam1", "diam2", "diam3","diam4","dist1",
                 "dist2", "dist3", "dist4" ,"species1", "species2", "species3", "species4", "az1", "az2", 
                 "az3", "az4", "corner", "sectioncorner","surveyyear", "point")]




# write as csv:
write.csv(used.data.alb, "data/outputs/Point_Data_From_Goringetal16_used_data_alb.csv", row.names = FALSE)

                       

