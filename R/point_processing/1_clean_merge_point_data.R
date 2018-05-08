## first step is to clean all the data for Southern MI, Uppermidwest, and Indiana + Illinois separately, get correction factors for all the data, then join together
## will work on combining the all the data before estimating correction factors, to make the correction factor generation more intuitive
## merge point data from UMW, IL, IN, SO MI

## TODO: add assertions at various places

## based on KAH 01a_clean_merge_IN_IL.r and Simon's [Jody/Kelly, what is file name from which this code was obtained]

library(readr)
library(dplyr)
library(fields)


final_columns <- c("x","y","twp","surveyyear",
                     "L1_tree1", "L1_tree2", "L1_tree3", "L1_tree4",
                     "L3_tree1", "L3_tree2", "L3_tree3", "L3_tree4",
                     "az1", "az2", "az3", "az4",
                     "dist1", "dist2", "dist3", "dist4",
                     "diam1", "diam2", "diam3", "diam4",
                     "corner", "sectioncorner","state", "point_id", "vegtype")

# ----------------------------------DATA CLEANING: IN + IL --------------------------------------------------

ind <- read_csv(file.path(raw_data_dir, indiana_file), guess_max = 100000)
il <- read_csv(file.path(raw_data_dir, illinois_file), guess_max = 100000)

ind <- ind %>% mutate(point_id = seq_len(nrow(ind)), vegtype = NA)
il <- il %>% mutate(point_id = seq_len(nrow(il)), vegtype = NA)

## Remove the following corners from IL because they are on the IL-WI border and these IL corners are very close to adjacent WI corners.  
## These IL corners are less than 165 m away from the WI corners and both sets of corners have similar taxa in that most are Oaks at Level 3a

IL_WI_overlap_points <- c("642067","642022","642315","642443","641958","642391","639084","639089","639087","639156","638543","639082","698794","696675","698842","699244","698820","696658")

il <- il %>% filter(!entry_id %in% IL_WI_overlap_points)

if(sum(is.na(ind$L1_tree1)) || sum(is.na(il$L1_tree1)))
    cat("Missing values in taxon for first tree in IN or IL.\n")

if(sum(is.na(ind$x)) || sum(is.na(ind$y) || sum(is.na(il$x)) || sum(is.na(il$y))))
    cat("Missing locations for points in IN or IL.\n")

if(sum(ind$L1_tree2 == "No data", na.rm = TRUE) || sum(il$L1_tree2 == "No data", na.rm = TRUE))
    cat("'No data' found for second tree in IN or IL; this case not handled by the code.\n")

## Our density/biomass calculations are on a per-land-area basis, excluding water area
ind <- ind %>% filter(!(L1_tree1 %in% c('No data', 'Water', 'Wet')))
il <- il %>% filter(!(L1_tree1 %in% c('No data', 'Water', 'Wet')))

## Converting L1 (survey abbreviation) to L3 (Paleon nomenclature) taxa; currently this overwrites existing L3 in Indiana (for Illinois it has already been removed) but ensuring use of current taxon conversion file; in future L3 will not be in the input files and will solely be created here.

## Indiana conversion
spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == indiana_conversion_domain) %>%
    select(level1, level3a)

ind <- ind %>% select(-L3_tree1, -L3_tree2, -L3_tree3, -L3_tree4) %>%  ## so we can replace these columns
    left_join(spec_codes, by = c('L1_tree1' = 'level1')) %>% rename(L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree2' = 'level1')) %>% rename(L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree3' = 'level1')) %>% rename(L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree4' = 'level1')) %>% rename(L3_tree4 = level3a)

if(sum(is.na(ind$L3_tree1)) != sum(is.na(ind$L1_tree1)) ||
   sum(is.na(ind$L3_tree2)) != sum(is.na(ind$L1_tree2)) ||
   sum(is.na(ind$L3_tree3)) != sum(is.na(ind$L1_tree3)) ||
   sum(is.na(ind$L3_tree4)) != sum(is.na(ind$L1_tree4)))
    cat("Apparently some L1 taxa are missing from the Indiana L1 to L3 conversion table.\n")

## Illinois conversion
spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == illinois_conversion_domain) %>%
    select(level1, level3a) 

il <- il %>% 
    left_join(spec_codes, by = c('L1_tree1' = 'level1')) %>% rename(L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree2' = 'level1')) %>% rename(L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree3' = 'level1')) %>% rename(L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree4' = 'level1')) %>% rename(L3_tree4 = level3a)

if(sum(is.na(il$L3_tree1)) != sum(is.na(il$L1_tree1)) ||
   sum(is.na(il$L3_tree2)) != sum(is.na(il$L1_tree2)) ||
   sum(is.na(il$L3_tree3)) != sum(is.na(il$L1_tree3)) ||
   sum(is.na(il$L3_tree4)) != sum(is.na(il$L1_tree4)))
    cat("Apparently some L1 taxa are missing from the Illinois L1 to L3 conversion table.\n")

## Note that check for NAs appearing in L3 that did not occur in L1 replaces a former
## replacement of NAs with 'No tree', which is not accurate since
## a non-available conversion is not the same as no tree appearing.
## If we really don't know the L3 taxon, we should convert to 'Unknown tree'.

## make sure township names have the state in front of them for later merging of state data:
ind <- ind %>% mutate(twp = paste0('IN_', TRP)) %>% dplyr::select(-TRP)
il <- il %>% mutate(twp = paste0('IL_', TRP)) %>% dplyr::select(-TRP)

if(any(ind$state != 'IN') || any(il$state != 'IL'))
    cat("State field missing from one or more rows in IN or IL.\n")

ind <- ind %>% rename(dist1 = chainstree, dist2 = chainstree2, dist3 = chainstree3, dist4 = chainstree4,
                      diameter1 = diameter, bearing1 = bearing, bearingdir1 = bearingdir, degrees1 = degrees)
il <- il %>% rename(dist1 = chainstree, dist2 = chainstree2, dist3 = chainstree3, dist4 = chainstree4,
                    diameter1 = diameter, bearing1 = bearing, bearingdir1 = bearingdir, degrees1 = degrees)

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
                     "cornerid", "typecorner","state", "point_id", "vegtype")

ind <- ind[columns_to_keep] 
il <- il[columns_to_keep]

inil <- rbind(ind, il)

inil <- inil %>% rename(diam1 = diameter1, diam2 = diameter2, diam3 = diameter3, diam4 = diameter4)

## Change 88888/99999 to NA but not in in bearing columns because an NA plus a degrees of 0 means a cardinal direction while 88888/9999 is unknown
## 88888/99999 in taxa have been dealt with in L1->L3 conversion
cols <- c("degrees1", "degrees2", "degrees3","degrees4",
          "dist1", "dist2", "dist3", "dist4",
          "diam1", "diam2", "diam3", "diam4")

inil[ , cols] <- sapply(inil[ , cols], convert_to_NA, missingCodes = c(88888,99999))


notree <- inil %>% filter(L1_tree1 == 'No tree')
if(sum(is.na(notree$dist1) & is.na(notree$dist2) & is.na(notree$dist3) & is.na(notree$dist4) &
       is.na(notree$diam1) & is.na(notree$diam2) & is.na(notree$diam3) & is.na(notree$diam4))
       != nrow(notree))
    cat("Found non-NA distances or diameters for no tree points in IN or IL.\n")

## create a survey year variable that corresponds to survey year correction factors

## We have some corners in IN & IL that are missing years. Based on exploratory analysis of their locations
## it is safe to assume these points were surveyed at a similar time as the points around them.
## This requires a bit over 1 GB RAM.

distances <- rdist(inil[inil$year == 9999, c('x', 'y')],
              inil[inil$year != 9999, c('x', 'y')])
closest <- apply(distances, 1, which.min)
inil$year[inil$year == 9999] <- inil$year[inil$year != 9999][closest]
if(sum(is.na(inil$year)) || min(inil$year < 1799) || max(inil$year > 1849))
    cat("Unexpected missing year or year outside 1799-1849 range.\n")
rm(distances)

## create a survey year variable that coresponds to survey year correction factors
## for IL/IN, '1825+' vs. '< 1825' and state uniquely defines correction factor - don't need region info 
inil <- inil %>% mutate(surveyyear = ifelse(year >= 1825, '>=1825', '<=1824'))

inil[ , paste0('az', 1:4)] <- get_angle_inil(as.matrix(inil[ , paste0('bearing', 1:4)]),
                           as.matrix(inil[ , paste0('degrees', 1:4)]))

if(max(inil[ , paste0('az', 1:4)], na.rm = TRUE) >= 360 | max(inil[ , paste0('az', 1:4)], na.rm = TRUE) < 0)
    cat("Found azimuths outside of 0-359 in IN/IL")

#----------Getting correction factors----------------------

##  Indiana and Illinois data have same correction factors for the whole state
##  Correction factors vary depending on which type of corner you are at
                       
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

inil <- inil %>% mutate(sectioncorner = ifelse(inil$cornerid %in% intsec | inil$cornerid %in% extsec, 
                                               'section', 'quartersection'), 
                        corner = ifelse(inil$cornerid %in% intsec | inil$cornerid %in% intqtr, 
                                        'internal', 'external'))

inil <- inil[final_columns]

# ----------------------------------DATA CLEANING: SOUTHERN MI --------------------------------------------------

somi <- read_csv(file.path(raw_data_dir, southern_michigan_file), guess_max = 100000)
somi <- somi %>% mutate(point_id = seq_len(nrow(somi)), vegtype = NA, state = "SoMI")

if(F) {
## Schoolcraft County and Isle Royale data provided separately but in same format as so MI data
nomi_extra <- read_csv(file.path(raw_data_dir, northern_michigan_supp_file), guess_max = 100000)
nomi_extra <- nomi_supp %>% mutate(point_id = seq_len(nrow(nomi_extra)), vegtype = NA, state = "NoMI_extra")
}

## These codes are not used in southern Michigan so don't need to do this filtering:
## somi <- somi %>% filter(!(species1 %in% c('No data', 'Water', 'Wet')))

## formerly we check for species1 and species2 being missing but having a positive diameter but
## there is only one case of missing species and existing diameter and that is 4th tree in
## case where it is not among the closest two. 

##  converting level 1 species to level 3 species:

## actual so MI points
spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == southern_michigan_conversion_domain) %>%
    select(level1, level3a) 

somi_lower <- somi %>% 
    left_join(spec_codes, by = c('species1' = 'level1')) %>% rename(L1_tree1 = species1, L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('species2' = 'level1')) %>% rename(L1_tree2 = species2, L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('species3' = 'level1')) %>% rename(L1_tree3 = species3, L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('species4' = 'level1')) %>% rename(L1_tree4 = species4, L3_tree4 = level3a)
if(F){ 
## Schoolcraft county and Isle Royale points
spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == upper_midwest_conversion_domain) %>%
    select(level1, level3a) 

nomi_extra <- nomi_extra %>% 
    left_join(spec_codes, by = c('species1' = 'level1')) %>% rename(L1_tree1 = species1, L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('species2' = 'level1')) %>% rename(L1_tree2 = species2, L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('species3' = 'level1')) %>% rename(L1_tree3 = species3, L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('species4' = 'level1')) %>% rename(L1_tree4 = species4, L3_tree4 = level3a)

somi <- rbind(somi_lower, nomi_extra)
} else {
    print("not yet reading Schoolcraft/IR")
    somi <- somi_lower 
}

if(sum(is.na(somi$L3_tree1)) != sum(is.na(somi$L1_tree1)) ||
   sum(is.na(somi$L3_tree2)) != sum(is.na(somi$L1_tree2)) ||
   sum(is.na(somi$L3_tree3)) != sum(is.na(somi$L1_tree3)) ||
   sum(is.na(somi$L3_tree4)) != sum(is.na(somi$L1_tree4)))
    cat("Apparently some L1 taxa are missing from the Michigan L1 to L3 conversion table.\n")

#remove the entries with sec_corner = "Check"
somi <- somi %>% filter(!sec_corner %in% 'Check')

## make sure township names have the state in front of them:
somi <- somi %>% mutate(twp = paste0('MI_', twnrng)) %>% mutate(state = 'MI')

## determine subdomain for use with correction factors
if(nrow(somi) != length(grep('[EW]', somi$range)))
    cat("Can't assign surveyyear for some southern Michigan sites.\n")
## shorthand (and unique) for 'SE - E of central Meridian S of tension'
surveyyear <- rep('<=1824', nrow(somi))  
## shorthand (and unique) for 'SW - W of central Meridian S of tension'
surveyyear[grep('W', somi$range)] <- '1825 -1835'  
## Schoolcraft and Isle Royale in UP: >=1840+ is shorthand (and unique) for 'UP, >=1840'
surveyyear[somi$point_y > 900000] <- '>=1840'

somi$surveyyear <- surveyyear

## There is some transcription error that shifted values between fields causing
## apparent large trees in SE MI. Suggested approach from Charlie Cogbill
## is to omit trees in SE MI with az1=0 or with non-missing 2nd tree and az2=0
## Note that in some cases this filters out points where the first two trees (by distance)
## have non-0 az values, but if we want to do based on two nearest, we would
## want to see if we are confident about distance values.
somi <- somi %>% filter(!(surveyyear == '<=1824' & az1 == 0 & !is.na(L1_tree1))) %>%
    filter(!(surveyyear == '<=1824' & az2 == 0 & !is.na(L1_tree2))) %>%
    filter(!(surveyyear == '<=1824' & az3 == 0 & !is.na(L1_tree3))) %>%
    filter(!(surveyyear == '<=1824' & az4 == 0 & !is.na(L1_tree4)))

## The az_360 columns were calculated using the quadrant and the az_x (x = 1:4), values.
## e.g., if the quadrant number was 1, the AZ as read from the mylar maps was used as is.
## If the quadrant number was 2, az1_360 is 180-az1.
## Quadrants of 3 were 180+az1 and quadrants of 4 were 360-az1.
somi <- somi %>% select(-az1, -az2, -az3, -az4) %>% 
    rename(az1 = az1_360, az2 = az2_360, az3 = az3_360, az4 = az4_360)  

somi <- somi %>% mutate(corner = ifelse(sec_corner == "Extsec", 'external', 'internal'),
                    sectioncorner = ifelse(cornertype == 'section',  'section', 'quartersection'))

somi <- somi %>% rename(x = point_x, y = point_y)

somi <- somi[final_columns]



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

if(!file.exists(file.path(raw_data_dir, minnesota_file))) 
    unzip(file.path(raw_data_dir, minnesota_zipfile), exdir = raw_data_dir)
if(!file.exists(file.path(raw_data_dir, wisconsin_file))) 
    unzip(file.path(raw_data_dir, wisconsin_zipfile), exdir = raw_data_dir)
if(!file.exists(file.path(raw_data_dir, northern_michigan_file))) 
    unzip(file.path(raw_data_dir, northern_michigan_zipfile), exdir = raw_data_dir)

wi <- read_csv(file.path(raw_data_dir, wisconsin_file), guess_max = 100000) %>% mutate(state = 'WI')
mn <- read_csv(file.path(raw_data_dir, minnesota_file), guess_max = 100000) %>% mutate(state = 'MN')
nomi <- read_csv(file.path(raw_data_dir, northern_michigan_file), guess_max = 100000) %>% mutate(state = 'NoMI')

wi <- wi %>% mutate(point_id = seq_len(nrow(wi)))
nomi <- nomi %>% mutate(point_id = seq_len(nrow(nomi)), vegtype = NA)
mn <- mn %>% mutate(point_id = seq_len(nrow(mn)))

names(wi) <- tolower(names(wi))
names(mn) <- tolower(names(mn))
names(nomi) <- tolower(names(nomi))

wi <- wi %>% rename(year = year_)

nomi <- nomi %>% rename(diam1 = dbh1, diam2 = dbh2, diam3 = dbh3, diam4 = dbh4,
                        sp1 = spp1, sp2 = spp2, sp3 = spp3, sp4 = spp4,
                        az1 = azimuth, az2 = azimuth2, az3 = azimuth3, az4 = azimuth4)

if(any(!wi$rangdir %in% c(0, 2, 4)))
    stop("Unexpected 'rangdir' values found in Wisconsin")

##  The wisconsin Range is set as a single value, the 'E' and 'W' codes are in
##  RANGDIR.  Looking at the data it also looks like there are a few ranges
##  that are miscoded with slight differences from cases Simon was checking. 
wi <- wi %>% filter(rangdir != 0) %>%   ## Single point that seems to have no data associated with it, so remove it.
    mutate(rangdir = ifelse(rangdir == '2', 'W', 'E')) %>%
    mutate(rangdir = ifelse(rangdir == 'W' & y_alb > 1050000 & x_alb > 4.6e5, 'E', rangdir)) %>%
    mutate(rangdir = ifelse(rangdir == 'E' & y_alb > 1050000 & x_alb < 4.3e5, 'W', rangdir)) %>%
    mutate(rng = paste0(as.character(range), rangdir)) %>%
    select(-range, -rangdir)

##  Michigan's point numbers are wrong in the dataset.  Simon is not sure where the 
##  error arose from, but we need them to be able to assign section & quartersection
##  points.
## Per discussion in github issue #31, the data come from two sources, so just
## need to use value from whichever is non-empty.
nomi$recnum = format(mi$recnum, scientific=F)
nomi$recnum_c = format(mi$recnum_c, scientific=F)
nomi <- nomi %>% mutate(newrecnum = pmax(recnum,recnum_c))
nomi <- nomi %>% mutate(point = as.numeric(substr(newrecnum, 9, 11)))

mn <- mn %>% rename(point = tic)


##  We have made a choice to say that all taxa labelled 'Beech' in Minnesota are likely
##  Bluebeech, or, in our dataset, Ironwood.
## Ideally this would be done in taxon translation table, but we don't have separate
## MN/WI/NoMI conversions
mn <- mn %>% mutate(sp1 = ifelse(sp1 == 'BE', 'IR', sp1),
                    sp2 = ifelse(sp2 == 'BE', 'IR', sp2),
                    sp3 = ifelse(sp3 == 'BE', 'IR', sp3),
                    sp4 = ifelse(sp4 == 'BE', 'IR', sp4))


##  We want the Minnesota data to reflect water in the same way that the Wisconsin data does.
##  Almendinger references the following land cover codes that are likely to have water:
##  'A' - Creek (unlikely to be exclusively water)
##  'M' - Marsh
##  'S' - Swamp (Charlie notes this meant wetland and not necessarily forested)
##  'L' - Lake
##  'R' - River
waterTypes <- c('L', 'M', 'S', 'R', 'A')

if(sum(mn$vegtype %in% waterTypes & mn$sp1 == '_' & (mn$sp2 != '_' | mn$sp3 != '_' | mn$sp4 != '_')))
    cat("MN water corners have trees with species info.\n")

mn <- mn %>% mutate(sp1 = ifelse(vegtype %in% waterTypes & sp1 == '_', 'QQ', sp1),
                    sp2 = ifelse(vegtype %in% waterTypes & sp2 == '_', 'QQ', sp2),
                    sp3 = ifelse(vegtype %in% waterTypes & sp3 == '_', 'QQ', sp3),
                    sp4 = ifelse(vegtype %in% waterTypes & sp4 == '_', 'QQ', sp4))

## a few of these cases have non-missing dist&diam; should check on these
mn <- mn %>% mutate(sp1 = convert_to_NA(sp1, '_'),
              sp2 = convert_to_NA(sp2, '_'),
              sp3 = convert_to_NA(sp3, '_'),
              sp4 = convert_to_NA(sp4, '_'))

## exclude water points -- all water as well as points with 1 tree, per issue #35
numQQ <- apply(mn[ , paste0('sp', 1:4)], 1, function(x) sum(x == 'QQ', na.rm = TRUE))
mn <- mn %>% filter(numQQ <= 2) 

## Next exclude no-tree points with unclear vegtype values

## these points occur mostly in northern Minnesota (Boundary Waters) as well as in
## Many of them are in east-west straight lines, suggesting not usable
mn <- mn %>% filter(!(is.na(sp1) & is.na(sp2) & is.na(sp3) & is.na(sp4) & vegtype == '_'))

## forest, grove, pine grove seem inconsistent with lack of trees, so exclude these points
forestedTypes <- c('F', 'G', 'J')
mn <- mn %>% filter(!(is.na(sp1) & is.na(sp2) & is.na(sp3) & is.na(sp4) & vegtype %in% forestedTypes))


mn_survey <- read_csv(file.path(raw_data_dir, minnesota_survey_file), guess_max = 10000)
mn_survey <- mn_survey %>% mutate(TOWN = paste0('T', formatC(TOWN, width=3, flag='0'), 'N'),
                                  RANG = paste0('R', formatC(RANG, width=2, flag='0'), RDIR)) %>%
    dplyr::select(TOWN, RANG, YEAR)
names(mn_survey) <- tolower(names(mn_survey))

mn <- mn %>% left_join(mn_survey, by = c('twp' = 'town', 'rng' = 'rang'))

## two points with twp/rng not in survey file
distances <- rdist(mn[is.na(mn$year), c('x_alb','y_alb')],
                   mn[!is.na(mn$year), c('x_alb','y_alb')])
closest <- apply(distances, 1, which.min)
mn$year[is.na(mn$year)] <- mn$year[!is.na(mn$year)][closest]

if(sum(is.na(mn$year)) || min(mn$year < 1847) || max(mn$year > 1907))
    cat("Unexpected missing year or year outside 1847-1907 range")

## '<=1853' vs. '>=1854' plus state uniquely defines correction factors without need for region info
mn <- mn %>% mutate(surveyyear = ifelse(year <= 1853, "<=1853", ">=1854"))


distances <- rdist(wi[wi$year == 0, c('x_alb','y_alb')],
                   wi[wi$year != 0, c('x_alb','y_alb')])
closest <- apply(distances, 1, which.min)
wi$year[wi$year == 0] <- wi$year[wi$year != 0][closest]

if(sum(is.na(wi$year)) || min(wi$year < 1832) || max(wi$year > 1891))
    cat("Unexpected missing year or year outside 1847-1907 range")

## '<=1845' vs. '>=1846' plus state uniquely defines correction factors without need for region info
wi <- wi %>% mutate(surveyyear = ifelse(year <= 1845, "<=1845", ">=1846"))

## 4088 cases; in essentially all of them, no info on trees 2-4, so assume these are fully water
wi <- wi %>% filter(sp1 != 'QQ')

## based on parsing notes field, do include some points that seem to be no-tree points
## (non-water and non-"tree is point/corner/post")
miss <- nomi %>% filter(is.na(sp1) & is.na(sp2) & is.na(sp3) & is.na(sp4))
nomi <- nomi %>% filter(!(is.na(sp1) & is.na(sp2) & is.na(sp3) & is.na(sp4)))

## points without any taxa and no notes are ambiguous so remove
miss <- miss %>% filter(!is.na(notes))
## water points
waterText <- c("(MARSH|POND|LAKE|WATER|RIVER|SWAMP|BROOK|STREAM|LK MICHIGAN)")
miss <- miss %>% filter(!grepl(waterText, notes, ignore.case = TRUE))
## points where tree is the corner so ambiguous what density would be
## some of these might be trees at previous points, but can't determine which
treeAsPostText <- c("(IS POST|IS CORNER|IS QUARTER|IS SECTION|AS CORNER|UPON CORNER|AS 1/4|IS 1/4|FOR 1/4|FOR CORNER|IS PSOT|ISPOST|AS POST|IS A QUARTER|IN CORNER|IS THE QUARTER CORNER|IS A QAURTER CORNER|IS A SECTION|A QUARTER CORNER|IS THE SECTION CORNER|FOR QUARTER|IS THE QAURTER CORNER)")
miss <- miss %>% filter(!grepl(treeAsPostText, notes, ignore.case = TRUE))
## indications that data lost or not noted
missingDataText <- c("(THEN LOST|INFORMATION NOT GIVEN|NOT IN NOTES|NO INFORMATION|RANDOM NOTES|OMITTED|OMMITTED|CUT OFF|MICROFILM|TREE IS DEAD|TA, RP, PS, WP|PS, RP|TRAIL COURSE)")
miss <- miss %>% filter(!grepl(missingDataText, notes, ignore.case = TRUE))

## otherwise, notes generally say 'no witness trees', 'no trees convenient', 'no bearing trees', 'no other tree data', 'no trees'; assumed to indicate no-tree points
nomi <- rbind(nomi, miss)


if(FALSE) { ## need to enable this once have code to split UP from northern LP
    ## >=1840 is shorthand (and unique) for 'UP, >=1840'
    ## >=1836 is shorthand (and unique) for 'north Lower, >=1836'
    nomi <- nomi %>% mutate(surveyyear = ifelse(domain == 'UP', '>=1840', '>=1836'))
}

nomi <- nomi %>% mutate(surveyyear = "allN")

wi <- wi %>% mutate(twp = paste0('WI_', township)) %>% select(-township)
mn <- mn %>% mutate(twp = paste0('MN_', twp))
nomi <- nomi %>% mutate(twp = paste0('MI_', twp))

columns_to_keep <- c("point", "twp", "rng", "surveyyear",
                     "sp1", "sp2", "sp3", "sp4",
                     "az1", "az2", "az3", "az4",
                     "dist1", "dist2", "dist3", "dist4",
                     "diam1", "diam2", "diam3", "diam4",
                     'x_alb', 'y_alb', 'state', 'point_id', 'vegtype')

mn <- mn[columns_to_keep] 
nomi <- nomi[columns_to_keep]
wi <- wi[columns_to_keep]

umw <- rbind(mn, wi, nomi)

umw <- umw %>% rename(x = x_alb, y = y_alb)


missingCodes <- c(8888, 9999)
azMissingCodes <- c('8888','9999','_')
umw <- umw %>% mutate(dist1 = convert_to_NA(dist1, missingCodes),
                      dist2 = convert_to_NA(dist2, missingCodes),
                      dist3 = convert_to_NA(dist3, missingCodes),
                      dist4 = convert_to_NA(dist4, missingCodes),
                      diam1 = convert_to_NA(diam1, missingCodes),
                      diam2 = convert_to_NA(diam2, missingCodes),
                      diam3 = convert_to_NA(diam3, missingCodes),
                      diam4 = convert_to_NA(diam4, missingCodes),
                      az1 = convert_to_NA(az1, azMissingCodes),
                      az2 = convert_to_NA(az2, azMissingCodes),
                      az3 = convert_to_NA(az3, azMissingCodes),
                      az4 = convert_to_NA(az4, azMissingCodes),
                      sp1 = convert_to_NA(sp1, 0),
                      sp2 = convert_to_NA(sp2, 0),
                      sp3 = convert_to_NA(sp3, 0),
                      sp4 = convert_to_NA(sp4, 0))

   

##  get_angle_um converts the four character azimuth (e.g. N43E) to a numeric, 360
##  degree angle.  It also has to deal with a number of special cases.
##  1186 cases left as NA after this processing - cases like 'N85', 'W45E'
umw[ , paste0('az', 1:4)] <- get_angle_umw(as.matrix(umw[ , paste0('az', 1:4)]))

spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
    filter(domain == upper_midwest_conversion_domain) %>%
    select(level1, level3a) 

umw <- umw %>% 
    left_join(spec_codes, by = c('sp1' = 'level1')) %>% rename(L1_tree1 = sp1, L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('sp2' = 'level1')) %>% rename(L1_tree2 = sp2, L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('sp3' = 'level1')) %>% rename(L1_tree3 = sp3, L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('sp4' = 'level1')) %>% rename(L1_tree4 = sp4, L3_tree4 = level3a)

if(sum(is.na(umw$L3_tree1)) != sum(is.na(umw$L1_tree1)) ||
   sum(is.na(umw$L3_tree2)) != sum(is.na(umw$L1_tree2)) ||
   sum(is.na(umw$L3_tree3)) != sum(is.na(umw$L1_tree3)) ||
   sum(is.na(umw$L3_tree4)) != sum(is.na(umw$L1_tree4)))
    cat("Apparently some L1 taxa are missing from the UMW L1 to L3 conversion table.\n")

#  Points within a township are either sections or quartersections.  This
#  is the list of points that are sections.  All others are quarter-sections.
sections <- c(2, 5, 8, 11, 14, 18, 21, 24, 27, 30,
              34, 37, 40, 43, 46, 50, 53, 56, 59, 62,
              66, 70, 74, 78, 82,
              87, 89, 91, 93, 95, 98, 100, 102, 104, 106, 108,
              109, 111, 113, 115, 117, 119, 122, 123, 124, 125, 126)

#  These are the points on the outside of each township.
external <- c(109:120, 97:108, 86:96, 122:126)


umw <- umw %>% mutate(corner = ifelse(point %in% external, 'external', 'internal'),
                      sectioncorner = ifelse(point %in% sections, 'section', 'quartersection'))

umw <- umw[final_columns]

mw <- rbind(umw, inil, somi)

#--------------Reorder the tree number by distance to the point-----------------

#  At this point we need to make sure that the species are ordered by distance
#  so that trees one and two are actually the closest two trees.

## Missing trees occur in UMW and seem to be dead trees.

## Most cases of "Missing" are in three clumps SW of Green Bay and have no taxa for other
## three trees and 0 dist and diam; throw these points out as they were not surveyed
## (Menominee lands), though the middle clump does not seem to have the 'XC' code one would expect
mw <- mw %>% filter(!(mw$L3_tree1 == "Missing" & !is.na(mw$L3_tree1) & 
                      is.na(mw$L3_tree2) & is.na(mw$L3_tree3) & is.na(mw$L3_tree4) &
                      mw$dist1 == 0 & mw$dist2 == 0 & mw$dist3 == 0 & mw$dist4 == 0))


## Treat remaining missing as scattered dead trees but do not set to NA because have
## dist/diam in general and don't want to induce 1-tree points
mw <- mw %>% mutate(L3_tree1 = ifelse(L3_tree1 == "Missing", "Unknown tree", L3_tree1),
                    L3_tree2 = ifelse(L3_tree2 == "Missing", "Unknown tree", L3_tree2),
                    L3_tree3 = ifelse(L3_tree3 == "Missing", "Unknown tree", L3_tree3),
                    L3_tree4 = ifelse(L3_tree4 == "Missing", "Unknown tree", L3_tree4))


## set dists to NA when there is not a tree there (NA, water, no tree) so
## zeroes are not interpreted as 0 distance
nontree_codes <- c("Water", "No tree")
mw <- mw %>% mutate(dist1 = ifelse(is.na(L3_tree1) | L3_tree1 %in% nontree_codes, NA, dist1),
                    dist2 = ifelse(is.na(L3_tree2) | L3_tree2 %in% nontree_codes, NA, dist2),
                    dist3 = ifelse(is.na(L3_tree3) | L3_tree3 %in% nontree_codes, NA, dist3),
                    dist4 = ifelse(is.na(L3_tree4) | L3_tree4 %in% nontree_codes, NA, dist4))

## per issue #39 we are calculating density and biomass based on trees below veil line
## with correction to get density for trees above the veil line
## set small trees (below 8 inch veil line) dists to Inf so we find the bigger trees as closest two
if(FALSE){
mw <- mw %>% mutate(diam1 = ifelse(diam1 < diameter_cutoff_inches, NA, diam1), 
                    diam2 = ifelse(diam2 < diameter_cutoff_inches, NA, diam2),
                    diam3 = ifelse(diam3 < diameter_cutoff_inches, NA, diam3),
                    diam4 = ifelse(diam4 < diameter_cutoff_inches, NA, diam4))


mw <- mw %>% mutate(dist1 = ifelse(diam1 < diameter_cutoff_inches, Inf, dist1), 
                    dist2 = ifelse(diam2 < diameter_cutoff_inches, Inf, dist2),
                    dist3 = ifelse(diam3 < diameter_cutoff_inches, Inf, dist3),
                    dist4 = ifelse(diam4 < diameter_cutoff_inches, Inf, dist4))
}
## TODO: probably set L3 to NA if diam is NA

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
    reorder_col_blocks('diam', ords) %>% 
    reorder_col_blocks('az', ords) 

## determine Pair (points where only two trees surveyed) vs 2nQ (four trees surveyed)
## used for correction factors (see issue #42)
## it would be rare to only find two trees if looking for four,
## and if we have zero or one tree, we don't use correction factors anyway
tmp <- as.matrix(mw[ , paste0('L3_tree', 1:4)])
tmp[tmp %in% c('No tree', 'Water')] <- NA
ntree <- apply(tmp, 1,function(x) sum(!is.na(x)))
mw <- mw %>% mutate(point = ifelse(ntree > 2, '2nQ', 'Pair'))


save(mw, file = 'cleaned_point.Rda')

