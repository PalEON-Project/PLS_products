## first step is to clean all the data for Southern MI, Uppermidwest, and Indiana + Illinois separately, get correction factors for all the data, then join together
## will work on combining the all the data before estimating correction factors, to make the correction factor generation more intuitive
## merge point data from UMW, IL, IN, SO MI

## based on KAH 01a_clean_merge_IN_IL.r and Simon's [Jody/Kelly, what is file name from which this code was obtained]

library(readr)
library(dplyr)
library(fields)

## all needed?
library(sp)
library(spdep)
library(rgdal)

# ----------------------------------DATA CLEANING: IN + IL --------------------------------------------------

ind <- read_csv(file.path(raw_data_dir, indiana_file), guess_max = 100000)
il <- read_csv(file.path(raw_data_dir, illinois_file), guess_max = 100000)

if(sum(is.na(ind$L1_tree1)) || sum(is.na(il$L1_tree1)))
    cat("Missing values in taxon for first tree in IN or IL.\n")

if(sum(is.na(ind$x)) || sum(is.na(ind$y) || sum(is.na(il$x)) || sum(is.na(il$y))))
    cat("Missing locations for points in IN or IL.\n")

if(sum(ind$L1_tree2 == "No data", na.rm = TRUE) || sum(il$L1_tree2 == "No data", na.rm = TRUE))
    cat("'No data' found for second tree in IN or IL; this case not handled by the code.\n")

## Our density/biomass calculations are on a per-land-area basis, excluding water area
ind <- ind %>% filter(!(L1_tree1 %in% c('No data', 'Water', 'Wet')))
il <- il %>% filter(!(L1_tree1 %in% c('No data', 'Water', 'Wet')))

## Converting L1 (survey abbreviation) to L3 (Paleon nomenclature) taxa; currently this overwrites existing L3 but ensuring use of current taxon conversion file; in future L3 will not be in the input files and will solely be created here.

spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == indiana_conversion_domain) %>%
    select(level1, level3a) 

ind <- ind %>% select(-L3_tree1, -L3_tree2, -L3_tree3, -L3_tree4) %>%  ## so can directly overwrite these columns
    left_join(spec_codes, by = c('L1_tree1' = 'level1')) %>% rename(L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree2' = 'level1')) %>% rename(L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree3' = 'level1')) %>% rename(L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree4' = 'level1')) %>% rename(L3_tree4 = level3a)

if(sum(is.na(ind$L3_tree1)) != sum(is.na(ind$L1_tree1)) ||
   sum(is.na(ind$L3_tree2)) != sum(is.na(ind$L1_tree2)) ||
   sum(is.na(ind$L3_tree3)) != sum(is.na(ind$L1_tree3)) ||
   sum(is.na(ind$L3_tree4)) != sum(is.na(ind$L1_tree4)))
    cat("Apparently some L1 taxa are missing from the Indiana L1 to L3 conversion table.\n")

spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == illinois_conversion_domain) %>%
    select(level1, level3a) 

il <- il %>% select(-L3_tree1, -L3_tree2, -L3_tree3, -L3_tree4) %>%  ## so can directly overwrite these columns
    left_join(spec_codes, by = c('L1_tree1' = 'level1')) %>% rename(L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree2' = 'level1')) %>% rename(L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree3' = 'level1')) %>% rename(L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree4' = 'level1')) %>% rename(L3_tree4 = level3a)

if(sum(is.na(ind$L3_tree1)) != sum(is.na(ind$L1_tree1)) ||
   sum(is.na(ind$L3_tree2)) != sum(is.na(ind$L1_tree2)) ||
   sum(is.na(ind$L3_tree3)) != sum(is.na(ind$L1_tree3)) ||
   sum(is.na(ind$L3_tree4)) != sum(is.na(ind$L1_tree4)))
    cat("Apparently some L1 taxa are missing from the Illinois L1 to L3 conversion table.\n")

## Note that check for NAs appearing in L3 that did not occur in L1 replaces a former
## replacement of NAs with 'No tree', which is not accurate since
## a non-available conversion is not the same as no tree appearing.
## If we really don't know the L3 taxon, we should convert to 'Unknown tree'.

## make sure township names have the state in front of them:
ind <- ind %>% mutate(twp = paste0('IN_', TRP)) %>% dplyr::select(-TRP)
il <- il %>% mutate(twp = paste0('IL_', TRP)) %>% dplyr::select(-TRP)

if(any(ind$state != 'IN') || any(il$state != 'IL'))
    cat("State field missing from one or more rows in IN or IL.\n")

ind <- ind %>% rename(dist1 = chainstree, dist2 = chainstree2, dist3 = chainstree3, dist4 = chainstree4, diameter1 = diameter, bearing1 = bearing, bearingdir1 = bearingdir, degrees1 = degrees)
il <- il %>% rename(dist1 = chainstree, dist2 = chainstree2, dist3 = chainstree3, dist4 = chainstree4, diameter1 = diameter, bearing1 = bearing, bearingdir1 = bearingdir, degrees1 = degrees)

ind$bearing1 <- paste(ind$bearing1, ind$bearingdir1, sep = '_')
ind$bearing2 <- paste(ind$bearing2, ind$bearingdir2, sep = '_')
ind$bearing3 <- paste(ind$bearing3, ind$bearingdir3, sep = '_')
ind$bearing4 <- paste(ind$bearing4, ind$bearingdir4, sep = '_')
il$bearing1 <- paste(il$bearing1, il$bearingdir1, sep = '_')
il$bearing2 <- paste(il$bearing2, il$bearingdir2, sep = '_')
il$bearing3 <- paste(il$bearing3, il$bearingdir3, sep = '_')
il$bearing4 <- paste(il$bearing4, il$bearingdir4, sep = '_')

columns_to_keep <- c("x","y","twp","year",
                     "L1_tree1", "L1_tree2", "L1_tree3", "L1_tree4",
                     "L3_tree1", "L3_tree2", "L3_tree3", "L3_tree4",
                     "bearing1", "bearing2", "bearing3", "bearing4",
                     "degrees1", "degrees2", "degrees3","degrees4",
                     "dist1", "dist2", "dist3", "dist4",
                     "diameter1", "diameter2", "diameter3", "diameter4",
                     "cornerid", "typecorner","state")


ind <- ind[columns_to_keep] 
il <- il[columns_to_keep]

inil <- rbind(ind, il)

## Don't change 88888/99999 in bearing or taxon columns as they are needed for angle calculations and L1->L3 conversion


cols <- c("degrees1", "degrees2", "degrees3","degrees4",
          "dist1", "dist2", "dist3", "dist4",
          "diameter1", "diameter2", "diameter3", "diameter4")

inil[ , cols] <- sapply(inil[ , cols], convert_to_NA)


notree <- inil %>% filter(L1_tree1 == 'No tree')
num_notree <- nrow(notree)
if(sum(is.na(notree$dist1) & is.na(notree$dist2) & is.na(notree$dist3) &
       is.na(notree$diameter1) & is.na(notree$diameter2) & is.na(notree$diameter3)) != num_notree)
    cat("Found non-NA distances or diameters for no tree points in IN or IL.\n")

## create a survey year variable that coresponds to survey year correction factors

## We have some corners in IN & IL that are missing years,
## which means we have to either discard b/c we dont know which correction factors to use,
## or impute the surveyyear.
## JP, KH, & CP determined that for most of these missing years, it is safe to
## assume these points were surveyed at a similar time as the points around them.
## Here we impute the year of these missing points:

## this is a lot faster than using gDistance
distances <- rdist(inil[inil$year == 9999, c('x', 'y')],
              inil[inil$year != 9999, c('x', 'y')])
closest <- apply(distances, 1, which.min)

inil$year[inil$year == 9999] <- inil$year[inil$year != 9999][closest]
if(sum(is.na(inil$year)) || min(inil$year < 1799) || max(inil$year > 1849))
    cat("Unexpected missing year or year outside 1799-1849 range")

# create a survey year variable that coresponds to survey year correction factors
inil <- inil %>% mutate(surveyyear = ifelse(year >= 1825, '1825+', '< 1825'))

## TODO: check this once get clarity from KAH on the get_angle function issues
inil[ , paste0('az', 1:4)] <- get_angle_inil(inil[ , paste0('bearing', 1:4)],
                           inil[ , paste0('degrees', 1:4)],
                           inil[ , paste0('dist', 1:4)])

## check this and use 'az'
inil <- inil %>% mutate(azimuths = ifelse(azimuths > 360, NA, azimuths))
# azimuths[which(azimuths > 360)] <- NA


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

inil$corner <- ifelse(inil$cornerid %in% intsec, 'intsec',
          ifelse(inil$cornerid %in% intqtr, 'intqtr',
          ifelse(inil$cornerid %in% extsec, 'extsec',
          ifelse(inil$cornerid %in% extqtr,  'extqtr',
          ifelse(inil$typecorner == "(1/4) Section", "intqtr",
          ifelse(inil$typecorner == "Section", "intsec", "extsec"))))))

## we should rename either 'typecorner' or 'cornertype' as the naming is confusing

inil <- inil %>% mutate(cornertype = paste0(corner, state),
                        internal = ifelse(!cornerid %in% c("extsec", "extqtr", "external"), 'internal', 'external'),
                        section =  ifelse(typecorner %in% "Section", 'section', 'quarter-section'),
                        point = rep("P", nrow(inil)))

inil <- inil %>% select(-cornerid, -typecorner)


## TODO: things we probably move to after combining all regions:
##  - reorder columns based on increasing distance
##  - deal with question of 2 trees at distance 1
##  - deal with azimuth for zero-distance trees

# ----------------------------------DATA CLEANING: SOUTHERN MI --------------------------------------------------

somi <- read_csv(file.path(raw_data_dir, southern_michigan_file), guess_max = 100000)

somi <- somi %>% filter(!(species1 %in% c('No data', 'Water', 'Wet')))

## make sure township names have the state in front of them:
somi <- somi %>% mutate(twp = paste0('MI_', town)) %>% dplyr::select(-town) %>%
    mutate(state = 'MI')

## CJP: I don't see these 4 trees.
if(FALSE) {
    somi[(is.na(somi$species1) & somi$diam1 > 0) | (is.na(somi$species2) & somi$diam2>0),] <- rep(NA, ncol(somi))  #  removes four records with no identified trees, but identified diameters
}


##  converting level 1 species to level 3 species:

spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == southern_michigan_conversion_domain) %>%
    select(level1, level3a) 

somi <- somi %>% select(-L3_tree1, -L3_tree2, -L3_tree3, -L3_tree4) %>%  ## so can directly overwrite these columns
    left_join(spec_codes, by = c('species1' = 'level1')) %>% rename(L1_tree1 = species1, L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('species2' = 'level1')) %>% rename(L1_tree2 = species2, L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('species3' = 'level1')) %>% rename(L1_tree3 = species3, L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('species4' = 'level1')) %>% rename(L1_tree4 = species4, L3_tree4 = level3a)

if(sum(is.na(somi$L3_tree1)) != sum(is.na(somi$L1_tree1)) ||
   sum(is.na(somi$L3_tree2)) != sum(is.na(somi$L1_tree2)) ||
   sum(is.na(somi$L3_tree3)) != sum(is.na(somi$L1_tree3)) ||
   sum(is.na(somi$L3_tree4)) != sum(is.na(somi$L1_tree4)))
    cat("Apparently some L1 taxa are missing from the Michigan L1 to L3 conversion table.\n")


## The az_360 columns were calculated using the quadrant and the az_x (x = 1:4), values.
## e.g., if the quadrant number was 1, the AZ as read from the mylar maps was used as is.
## If the quadrant number was 2, az1_360 is 180-az1.
## Quadrants of 3 were 180+az1 and quadrants of 4 were 360-az1.
## NOTE: Missing azimuths are listed as "0" in az_360.

## in so. mi, all the azimuths == 0 are actually NA (from communication with Charlie):

somi <- somi %>% mutate(az1 = convert_to_NA(az1_360, 0),
                    az2 = convert_to_NA(az2_360, 0),
                    az3 = convert_to_NA(az3_360, 0),
                    az4 = convert_to_NA(az4_360, 0)) %>%
    select(-az1_360, -az2_360, -az3_360, -az4_360)
                    


## determine subdomain for use with correction factors
## should we find domain of closest point for those without 'E' or 'W'
surveyyear <- rep(NA, nrow(somi))
surveyyear[grep('E', somi$twnrng)] <- "SE"
surveyyear[grep('W', somi$twnrng)] <- "SW"
somi$surveyyear <- surveyyear

somi <- somi %>% rename(corner = sec_corner)

num_notrees <- (somi$species1 == 'No tree') + (somi$species2 == 'No tree')

somi <- somi %>% mutate(corner = ifelse(corner == "Extsec", 'external', 'internal'),
                    point = ifelse(num_notrees == 2, 'P', '2nQ'),
                    sectioncorner = ifelse(corner == 'section',  'section', 'quarter-section'))

somi <- somi %>% rename(x = point_x, y = point_y, twp = twnrng)

                       
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

unzip(file.path(raw_data_dir, upper_midwest_file), exdir = raw_data_dir)

## check on whether dist/diam are non-numeric

wisc <- readOGR(file.path(raw_data_dir, wisconsin_file), stringsAsFactors = FALSE)
minn <- readOGR(file.path(raw_data_dir, minnesota_file), stringsAsFactors = FALSE)
mich <- readOGR(file.path(raw_data_dir, michigan_file), stringsAsFactors = FALSE)

#  The files are in unique projections, this standardizes the projections to
#  a long-lat projection:
## Do we do anything with this?
if(FALSE) {
wisc <- spTransform(wisc, CRS('+proj=longlat +ellps=WGS84'))
minn <- spTransform(minn, CRS('+proj=longlat +ellps=WGS84'))
mich <- spTransform(mich, CRS('+proj=longlat +ellps=WGS84'))
}

clean_df <- function(input) {
    input <- as.data.frame(input)
    names(input) <- tolower(names(input))
    return(as_tibble(input))
}

wi <- clean_df(wisc)
mn <- clean_df(minn)
nomi <- clean_df(mich)

if(any(!wi$rangdir %in% c(0, 2, 4)))
    stop("Unexpected 'rangdir' values found in Wisconsin")

## recheck this next bit with WI CSV 

##  The wisconsin Range is set as a single value, the 'E' and 'W' codes are in
##  RANGDIR.  Looking at the data it also looks like there are a few ranges
##  that are miscoded (Simon had the W + >5e5 cases but none seen now, while
##  Simon did not have the W + < 4.8e5 cases that are seen now.
wi <- wi %>% filter(rangdir != 0) %>%   ## Single point that seems to have no data associated with it, so remove it.
    mutate(rangdir = ifelse(rangdir == '2', 'W', 'E')) %>%
    mutate(rangdir = ifelse(rangdir == 'W' & coords.x1 > 5e5, 'E', rangdir)) %>%
    mutate(rangdir = ifelse(rangdir == 'E' & coords.x2 < 4.8e5, 'W', rangdir)) %>%
    mutate(rng = paste0(as.character(range), rangdir)) %>%
    select(-range, -rangdir)

wi <- wi %>% mutate(twp = paste0('WI_', township)) %>% select(-township)
mn <- mn %>% mutate(twp = paste0('MN_', twp))
nomi <- nomi %>% mutate(twp = paste0('MI_', township))

##  Michigan's point numbers are wrong in the dataset.  Simon is not sure where the 
##  error arose from, but we need them to be able to assign section & quartersection
##  points.
nomi <- nomi %>% mutate(pnt = as.numeric(substr(recnum_c, 9, 11)))

mn <- mn %>% rename(point = tic)
nomi <- nomi %>% rename(point = pnt, sp1 = spp1, sp2 = spp2, sp3 = spp3, sp4 = spp4,
                        diam1 = dbh1, diam2 = dbh2, diam3 = dbh3, diam4 = dbh4,
                        az1 = azimuth, az2 = azimuth2, az3 = azimuth3, az4 = azimuth4)


## check on character/numeric issue
if(FALSE) {
    minn$DIST1 <- as.numeric(levels(minn$DIST1)[minn$DIST1])
    minn$DIST2 <- as.numeric(levels(minn$DIST2)[minn$DIST2])
    minn$DIST3 <- as.numeric(levels(minn$DIST3)[minn$DIST3])
    minn$DIST4 <- as.numeric(levels(minn$DIST4)[minn$DIST4])
}


##  We have made a choice to say that all taxa labelled 'Beech' in Minnesota are likely
##  Bluebeech, or, in our dataset, Ironwood.
## Ideally this would be done in taxon translation table, but we don't have separate
## MN/WI/NoMI conversions
mn <- mn %>% mutate(sp1 = ifelse(sp1 == 'BE', 'IR', sp1),
                    sp2 = ifelse(sp2 == 'BE', 'IR', sp2),
                    sp3 = ifelse(sp3 == 'BE', 'IR', sp3),
                    sp4 = ifelse(sp4 == 'BE', 'IR', sp4),


##  We want the Minnesota data to reflect water in the same way that the Wisconsin data does.
##  Almendinger references the following land cover codes that are likely to have water:
##  'A' - Creek (unlikely to be exclusively water)
##  'M' - Marsh
##  'S' - Swamp
##  'L' - Lake
##  'R' - not sure what this is - not indicated in Simon's code
waterTypes <- c('L', 'M', 'S', 'R', 'A')
mn <- mn %>% mutate(sp1 = ifelse(vegtype %in% waterTypes & sp1 == '_', 'QQ', sp1),
                    sp2 = ifelse(vegtype %in% waterTypes & sp1 == '_', 'QQ', sp1),
                    sp3 = ifelse(vegtype %in% waterTypes & sp1 == '_', 'QQ', sp1),
                    sp4 = ifelse(vegtype %in% waterTypes & sp1 == '_', 'QQ', sp1))

if(sum(mn$vegtype %in% waterTypes & mn$sp1 == '_' & (mn$sp2 != '_' | mn$sp3 != '_' | mn$sp4 != '_')))
    cat("MN water corners have trees with species info.\n")


mn_survey <- read_csv(file.path(raw_data_dir, minnesota_survey_file), guess_max = 10000)
mn_survey <- mn_survey %>% mutate(TOWN = paste0('T', formatC(TOWN, width=3, flag='0'), 'N'),
                                  RANG = paste0('R', formatC(RANG, width=2, flag='0'), RDIR)) %>%
    dplyr::select(TOWN, RANG, YEAR)

mn <- mn %>% left_join(mn_survey, by = c('twp' = 'TOWN', 'rng' = 'RANG'))
mn <- mn %>% rename(surveyyear = YEAR)

## two points with twp/rng not in survey file
distances <- rdist(mn[is.na(mn$surveyyear), c('x','y')],
                   mn[!is.na(mn$surveyyear), c('x','y')])
closest <- apply(distances, 1, which.min)

mn$surveyyear[is.na(mn$surveyyear)] <- mn$surveyyear[!is.na(mn$surveyyear)][closest]

if(sum(is.na(mn$surveyyear)) || min(mn$surveyyear < 1847) || max(mn$surveyyear > 1907))
    cat("Unexpected missing year or year outside 1847-1907 range")

mn <- mn %>% mutate(surveyyear = ifelse(surveyyear <= 1855, "1855", "1907"))



#  We need to bin the year information so that we can use it to calculate
#  appropriate Cottam Correction factors.  The survey instructions for the PLS
#  change at a number of points during the sirveys in Wisconsin, but are
#  considered to be fixed by the time.
#  Some wisconsin samples don't have a year.  Look this up and figure out why.
#  It causes a problem with the Cottam correction factor.

## the inequalities here are not precise, but this is what we have from Simon
## only current ambiguity is whether 1846 should be assigned to "1845"
wi <- wi %>% mutate(surveyyear = ifelse(year_ > 1851, '1851+',
                                 ifelse(year_ > 1846, '1846-1851',
                                 ifelse(year_ > 1834, '1834-1846',
                                 ifelse(year_ >= 1832, '1832-1834', NA)))))

distances <- rdist(wi[is.na(wi$surveyyear), c('x','y')],
                   wi[!is.na(wi$surveyyear), c('x','y')])
closest <- apply(distances, 1, which.min)

wi$surveyyear[wi$surveyyear)] <- wi$surveyyear[!is.na(wi$surveyyear)][closest]

if(sum(is.na(wi$surveyyear)) || min(wi$surveyyear < 1832) || max(wi$surveyyear > 1891))
    cat("Unexpected missing year or year outside 1847-1907 range")

wi <- wi %>% mutate(surveyyear = ifelse(surveyyear %in% c('1832-1834', '1834-1846'), "1845", "1907"))

#  There are also some weird Michigan points:
#  1.  Michigan has a set of points with NA as SPP1 but identifiable trees listed as
#      'tree'.  1549 of these are quartersection points, 45 are section points.  This
#      is clearly an artifact of the sampling method.  We remove these points.
#  2.  There are also 2909 no tree points in Michigan.  Most of these points are quarter
#      section points, and there is clear grographic bias.  We assume these points are
#      early survey points and remove them entirely.
## CJP: case 1 is a subset of case 2, so just remove case 2

nomi <- nomi %>% filter(!is.na(sp1))
## check this: there are ~10 cases where sp2,sp3, or sp4 not NA but sp1 is NA

if(FALSE) {  ## We computed these values, but then assigned all No MI data to one correction factor
    nomi <- nomi %>% mutate(surveyyear = ifelse(survyr > 1851, '1851+',
                                         ifelse(survyr > 1846, '1846-1851',
                                         ifelse(survyr > 1834, '1834-1846',
                                         ifelse(survyr >= 1832, '1832-1834', NA)))))
    
    distances <- rdist(nomi[is.na(nomi$surveyyear), c('x','y')],
                       nomi[!is.na(nomi$surveyyear), c('x','y')])
    closest <- apply(distances, 1, which.min)
    
    nomi$surveyyear[is.na(nomi$surveyyear)] <- nomi$surveyyear[!is.na(nomi$surveyyear)][closest]

    if(sum(is.na(nomi$surveyyear)) || min(nomi$surveyyear < 1836) || max(nomi$surveyyear > 1888))
        cat("Unexpected missing year or year outside 1836-1888 range")
}

nomi <- nomi %>% mutate(surveyyear = "allN")


# what about twp, rng etc?
columns_to_keep <- c("point", "twp", "rng", "surveyyear",
                     "sp1", "sp2", "sp3", "sp4",
                     "az1", "az2", "az3", "az4",
                     "dist1", "dist2", "dist3", "dist4",
                     "diam1", "diam2", "diam3", "diam4")


mn <- mn[columns_to_keep] 
nomi <- nomi[columns_to_keep]
wi <- wi[columns_to_keep]

umw <- rbind(mn, wi, nomi)

#  There are a set of 9999 values for distances which I assume are meant to be NAs.  Also, there are a set of points where
#  the distance to the tree is 1 or 2 feet.  They cause really big density estimates!
umw@data [ umw@data == '9999'] <- NA
umw@data [ umw@data == '8888'] <- NA
umw@data [ umw@data == '_'] <- NA       # Except those that have already been assigned to 'QQ'
umw@data [ umw@data == '99999'] <- NA
umw@data [ umw@data == '999999'] <- NA
umw@data [ umw@data == '6666'] <- NA
umw@data [ umw@data == '999'] <- NA
umw$DIAM1[is.na(umw$DIAM1)] <- 0
umw$DIAM2[is.na(umw$DIAM2)] <- 0




#  getAngle converts the four character azimuth (e.g. N43E) to a numeric, 360
#  degree angle.  It also has to deal with a number of special cases.
#  The code for getAngles is a bit scuzzy, but it leaves only 231 azimuths 
#  untranslated, this is a manageable number.
## TODO: check this once get clarity from KAH on the get_angle function issues
umw[ , paste0('az', 1:4)] <- get_angle(umw[ , paste0('az', 1:4)])

## azimuths <- apply(azimuths, 2, get_angle)


spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == upper_midwest_conversion_domain) %>%
    select(level1, level3a) 

umw <- umw %>% 
    left_join(spec_codes, by = c('sp1' = 'level1')) %>% rename(L1_tree1 = species1, L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('sp2' = 'level1')) %>% rename(L1_tree2 = species2, L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('sp3' = 'level1')) %>% rename(L1_tree3 = species3, L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('sp4' = 'level1')) %>% rename(L1_tree4 = species4, L3_tree4 = level3a)


#  Points within a township are either sections or quartersections.  This
#  is the list of points that are sections.  All others are quarter-sections.
sections <- c(2, 5, 8, 11, 14, 18, 21, 24, 27, 30,
              34, 37, 40, 43, 46, 50, 53, 56, 59, 62,
              66, 70, 74, 78, 82,
              87, 89, 91, 93, 95, 98, 100, 102, 104, 106, 108,
              109, 111, 113, 115, 117, 119, 122, 123, 124, 125, 126)

#  These are the points on the outside of each township.
external <- c(109:120, 97:108, 87, 89, 91, 93, 95, 122:126)


umw <- umw %>% filter(corner = ifelse(point %in% external, 'external', 'internal'),
                      sectioncorner = ifelse(point %in% sections, 'section', 'quartersection'))


mw <- rbind(umw, inil, somi)

#--------------Reorder the tree number by distance to the point-----------------

#  At this point we need to make sure that the species are ordered by distance
#  so that trees one and two are actually the closest two trees.

## find reordering vector for each row
ords <- t(apply(as.matrix(mw[ , paste0('dist', 1:4)]), 1, order, na.last = TRUE))

reorder_col_blocks <- function(data, colname, ords) {
    ## applies reordering vector by row to a 4-column block
    rowIndex <- 1:nrow(data)
    cols <- paste0(colname, 1:4)
    tmp <- as.matrix(data[ , cols])
    for(j in 1:4)
        data[ , cols[j]] <- tmp[cbind(rowIndex, ords[ , j])]
    return(data)
}

mw <- mw %>%
    reorder_col_blocks('dist', ords) %>% 
    reorder_col_blocks('L1_tree', ords) %>% 
    reorder_col_blocks('L3_tree', ords) %>% 
    reorder_col_blocks('bearing', ords) %>% 
    reorder_col_blocks('degrees', ords) %>%
    reorder_col_blocks('diameter', ords)


## remove 1-tree 0-distance points

## keep 2-tree points regardless of distances and truncate density (at say 1000 for now)

## keep 1-tree points and set density to something like 1 tree per unit circle of radius (150 links for now)

## set azimuth to 0 when distance is zero
## (This is being done for the moment because Simon originally
## needed a non-NA value in these cases; not clear if this is still needed

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

## CJP note: I don't see any '' cases; as far as NA, I think we
## want to leave as NA for now rather than treat as No tree
#  We need to indicate water and remove it.  There are 43495 cells with 'water'
#  indicated, and another 784 cells with 'missing' data.
#  when we limit these to the first two columns of the species table we get
#  a total of 25416 samples removed.

#  There are a set of dead taxa (DA, DB & cetera) that we exclude.  Only AM is
#  unknown at this point.  This excludes 213 trees.
species[species %in% ''] <- 'No tree'

#  Now we assign species that don't fit to the 'No tree' category.
species[is.na(species)] <- 'No tree'

# There is some cleaning to do.  A bit frustrating.  We can't confirm the diameters of
#  a number of points, although we hope to at some point in the future:
#  No stem density removals, none of the plots look like they have 'weird' points.
#  Basal area removals:
umw@data[which(as.numeric(umw$DIAM1) >100),] <- rep(NA, ncol(umw))  #  removes 19 trees with reported diameters over 250cm.
umw@data[which(as.numeric(umw$DIAM2) >100),] <- rep(NA, ncol(umw))  #  removes an additional 14 trees.
umw@data[(is.na(umw$SP1) & umw$DIAM1>0) | (is.na(umw$SP2) & umw$DIAM2>0),] <- rep(NA, ncol(umw))  #  removes four records with no identified trees, but identified diameters

#  When the object is NA, or the species is not a tree (NonTree or Water), set
#  the distance to NA.
dists[is.na(species) | species %in% c('No tree', 'Water', 'Missing')] <- NA


## check for 'Missing' as taxon code

# now kill missing cells:
final.data <- final.data[!final.data$level3a_1 %in% c('Water', 'Missing'),] 
final.data <- final.data[!final.data$level3a_2 %in% c('Water', 'Missing'),]
