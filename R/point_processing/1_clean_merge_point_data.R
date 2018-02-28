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
#  There are 8 corners that have 4 trees, and some corners with 3 trees. 
#  For now, we are treating all corners as '2 trees' for the correction factors
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

corner<- ifelse(inil$cornerid %in% intsec, 'intsec',
       ifelse(inil$cornerid %in% intqtr, 'intqtr',
              ifelse(inil$cornerid %in% extsec, 'extsec',
                     ifelse(inil$cornerid %in% extqtr,  'extqtr',
                            ifelse(inil$typecorner == "(1/4) Section", "intqtr",
                                   ifelse(inil$typecorner == "Section", "intsec", "extsec"))))))
       

inil$cornertype <- paste0(corner, inil$state)


#These are the columns for the final dataset.

final.data <- data.frame(inil$x,
                         inil$y,
                         inil$twp,
                         as.character(inil$state),
                         ranked.data[,1:8],
                         species[,1:4],
                         ranked.data[,13:16], 
                         inil$cornertype,
                         inil$surveyyear,
                         stringsAsFactors = FALSE)

colnames(final.data) <- c('PointX','PointY', 'Township','state',
                          paste('diam',    1:4, sep =''),
                          paste('dist',    1:4, sep = ''), 
                          paste('species', 1:4, sep = ''),
                          paste('az',      1:4, sep = ''), 'corner', 'surveyyear')

                          
Pair <- paste0(as.character(full.final$corner), full.final$surveyyear)

corr.vals <- read.csv('data/charlie_corrections_full_midwest.csv') # csv with correction factors from charlie--needs to be versioned properly

# get internal vs. ext and qtr vs section corners:                       
internal <- ifelse(!inil$cornerid %in% c("extsec", "extqtr"), 'internal', 'external')
#trees    <- ifelse(plot.trees == 2, 'P', '2NQ')
section  <- ifelse(inil$typecorner %in% "Section", 'section', 'quarter-section')
state <- final.data$state

# indiana and illinois are based on survey year:                       
corr.year     <- as.character(final.data$surveyyear)

# match up the correction facters with the inil corner types:
match.vec <- apply(corr.vals[,c("Pair", "year", "corner", "sectioncorner")], 1, paste, collapse = '')
to.match <- apply(data.frame(state, corr.year, internal, section, stringsAsFactors = FALSE), 1, paste, collapse = '')

corrections <- corr.vals[match(to.match, match.vec),]


#write the data & the density correction factors as a csv
write.csv(corrections, 'data/correction_factors.csv')
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

# make a dataframe with the values from Q1, Q2, Q3, Q4:
qvals <- cbind(as.numeric(nwmw$Q1), 
               as.numeric(nwmw$Q2),
               as.numeric(nwmw$Q3),
               as.numeric(nwmw$Q4))


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


#  These are the columns for the final dataset.

final.data <- data.frame(nwmw$point_x,
                         nwmw$point_y,
                        nwmw$twnrng,
                        state,
                        ranked.data[,1:8],
                        species,
                        ranked.data[,13:16],
                        corner,
                        nwmw$cornertype,
                        year ,
                        nwmw$twnrng,
                        stringsAsFactors = FALSE)

colnames(final.data) <- c('PointX','PointY', 'Township',"state",
                          paste('diam',    1:4, sep =''),
                          paste('dist',    1:4, sep = ''), 
                          paste('species', 1:4, sep = ''),
                          paste('az',      1:4, sep = ''), 'corner',"cornertype",'surveyyear')

# charlie designated all the NA azimuths as '0', and true 0 azimuths == 360, so get rid of the 0 values here:
final.data$az1[final.data$az1 <= 0 ] <- NA
final.data$az2[final.data$az2 <= 0] <- NA
final.data$az3[final.data$az3 <= 0] <- NA
final.data$az4[final.data$az4 <= 0] <- NA


final.data <- final.data[!is.na(final.data$PointX),]
# kill ths cells that are not == Extsec or ==Intsec

final.data <- final.data[final.data$corner %in% c("Extsec", "Intsec"),]

# now kill missing cells:

# there are a few strange points with erroneous X or Y values. Get rid of them here:
final.data <- final.data[ !final.data$PointX < 1000, ]
final.data <- final.data[ !final.data$PointY < 1000, ]


#write data to a csv:
write.csv(final.data, "data/lower_mi_final_data.csv")
#note there are still many NA values in the dataset--need to remove these!

used.data <- final.data
                       
# --------------------generate correction factors for southern mi------------------------------:
# read in the correction factors provided by Charlie Cogbill:
# correction factors for southern MI are based on the spatial location of the survey point + the corner type
corr.vals <- read.csv('data/charlie_corrections_full_midwest.csv')
correction <- data.frame(kappa = rep(NA, length(used.data)),
                         theta = rep(NA, length(used.data)),
                         zeta  = rep(NA, length(used.data)),
                         phi   = rep(NA, length(used.data)))

species2table <- data.frame(species1 = final.data$species1,
                            species2 = final.data$species2,
                            species3 = final.data$species3,
                            species4 = final.data$species4)
plot.trees <- rowSums(!(species2table == 'Water' | species2table == 'No tree'), na.rm = TRUE)

point.no <- as.character(final.data$corner)

#  So there are a set of classes here, we can match them all up:

internal <- ifelse(!point.no %in% "Extsec", 'internal', 'external')
trees    <- ifelse(plot.trees == 2, 'P', '2nQ')
section  <- ifelse(final.data$cornertype %in% "section", 'section', 'quarter-section')
state <- final.data$state


corr.year     <- as.character(final.data$year)
corr.year[state == 'MI' & final.data$Township %like% "W"] <- 'SW'
corr.year[state == 'MI' & final.data$Township %like% "E"] <- 'SE'

match.vec <- apply(corr.vals[,c("Pair", "year", "corner", "sectioncorner", "point")], 1, paste, collapse = '')
to.match <- apply(data.frame(state, corr.year, internal, section,trees, stringsAsFactors = FALSE), 1, paste, collapse = '')

correction <- corr.vals[match(to.match, match.vec),]


write.csv(correction, 'data/MI_correction_factors.csv')


                       
# ----------------------------------DATA CLEANING: UMW -------------------------------------------------------
 # data cleaning done in simon's witness tree code
 # correction factors for UMW are also generated in the witness tree code
  
                       
                       
                       
#-----------------------Merge all data and correction factors from UMW, IL, IN, SO MI------------------------------------------------

final.data <- read.csv(paste0("outputs/ndilin_pls_for_density_v",version,".csv"), stringsAsFactors = FALSE)
# corrections for stem density:
correction.factor <- read.csv("data//correction_factors.csv", header = TRUE)
colnames(correction.factor) <- c("X","Pair","regions","year","corner" ,      
                                "sectioncorner", "point" ,"Qdrt.model","kappa","theta" ,       
                                  "zeta","phi")
final.data <- final.data[,1:23]

# read in final data from michigan
final.data.mi <- read.csv(paste0("data/lower_mi_final_data.csv"), stringsAsFactors = FALSE)
final.data.mi <- final.data.mi[!names(final.data.mi) %in% c("cornertype", "NA.")]

# corrections for stem density:
correction.factor.mi <- read.csv("data//MI_correction_factors.csv", header = TRUE)


# read in final data from UMW 
# final.data.umw <- 
# corrections.umw <- 

# add the lower MI data below the INIL data: 

final.data <- rbind(final.data, final.data.mi, final.data.uwm)
correction.factor <- rbind(correction.factor, correction.factor.mi, correction.factor.umw)

