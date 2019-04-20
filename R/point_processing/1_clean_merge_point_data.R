## Clean all the data for southern MI, Upper Midwest (Minnesota, Wisconsin, northern Michigan),
## and {Indiana + Illinois + Detroit} separately,
## then combine and ensure that trees are ordered by distance.

## Run time for this file: approximately 1 minute

library(readr)
library(dplyr)
library(fields)
library(assertthat)

## common column names for merging data subsets
final_columns <- c("x","y","twp","surveyyear",
                     "L1_tree1", "L1_tree2", "L1_tree3", "L1_tree4",
                     "L3_tree1", "L3_tree2", "L3_tree3", "L3_tree4",
                     "az1", "az2", "az3", "az4",
                     "dist1", "dist2", "dist3", "dist4",
                     "diam1", "diam2", "diam3", "diam4",
                     "corner", "sectioncorner","state", "point_id", "vegtype")

## -------------------------DATA CLEANING: IN + IL + Detroit area ----------------------------------------------

## Detroit area was re-entered by Notre Dame so in same format as IL/IN

if(!file.exists(file.path(raw_data_dir, indiana_file))) 
    unzip(file.path(raw_data_dir, indiana_zipfile), exdir = raw_data_dir)
if(!file.exists(file.path(raw_data_dir, illinois_file))) 
    unzip(file.path(raw_data_dir, illinois_zipfile), exdir = raw_data_dir)
if(!file.exists(file.path(raw_data_dir, detroit_file))) 
    unzip(file.path(raw_data_dir, detroit_zipfile), exdir = raw_data_dir)

ind <- read_csv(file.path(raw_data_dir, indiana_file), guess_max = 100000)
il <- read_csv(file.path(raw_data_dir, illinois_file), guess_max = 100000)
det <- read_csv(file.path(raw_data_dir, detroit_file), guess_max = 100000)

ind <- ind %>% mutate(point_id = seq_len(nrow(ind)), vegtype = NA)
il <- il %>% mutate(point_id = seq_len(nrow(il)), vegtype = NA)
## 'state' needs to be 'Detroit' so don't have point_id values that overlap with point_id
## values for other parts of Michigan
det <- det %>% mutate(point_id = seq_len(nrow(det)), vegtype = NA, state = "Detroit")

## If all NAs, read in as logical and this causes problems with joining to spec_codes later.
if(class(det$L1_tree3) == "logical")
    det <- det %>% mutate(L1_tree3 = as.character(L1_tree3))
if(class(det$L1_tree4) == "logical")
    det <- det %>% mutate(L1_tree4 = as.character(L1_tree4))

## Remove the following corners from IN because they are 13 sets of two corners with identical tree information. 
## Jody has checked the survey notes for each set of corners and the tree information is identical. 
## We expect that the surveyors or the original transcribers in the General Land Office made errors in the transcription. 
## Since we do not know which corners the tree information comes from we are removing data from both corners

IN_duplicates <- c("662751","662755","R_1LJGyvH7k2ihanj","R_1HUnqzw2giTgXFr","600423","600396","82442","82457","713965","713968",
                   "704640","704623","705088","705091","81949","81951","605089","605085","R_cZ5Z28lG2ZhD7Ex","R_3Db6wD9TJdiFlCB",
                   "661407","661405","659856","659848","608889","603781")

ind <- ind %>% filter(!entry_id %in% IN_duplicates)

## Remove the following corners from IL because they are 2 sets of two corners with identical tree information. 
## Remove these 4 corners for the same reason as the IN duplicates above

IL_duplicates <- c("648480","648496","884","820")

il <- il %>% filter(!entry_id %in% IL_duplicates)

## Remove the following corners from IL because they are on the IL-WI border and these IL corners are very close to adjacent WI corners.  
## These IL corners are less than 165 m away from the WI corners and both sets of corners have similar taxa in that most are Oaks at Level 3a

IL_WI_overlap_points <- c("642067","642022","642315","642443","641958","642391","639084","639089","639087","639156","638543","639082","698794","696675","698842","699244","698820","696658")

il <- il %>% filter(!entry_id %in% IL_WI_overlap_points)

assert_that(!sum(is.na(ind$L1_tree1)) && !sum(is.na(il$L1_tree1)) && !sum(is.na(det$L1_tree1)),
            msg = "Missing values in taxon for first tree in IN or IL or Detroit.")

assert_that(!sum(is.na(ind$x)) && !sum(is.na(ind$y)) &&
            !sum(is.na(il$x)) && !sum(is.na(il$y)) &&
            !sum(is.na(det$x)) && !sum(is.na(det$y)),
    msg = "Missing locations for points in IN or IL or Detroit.")

nodata_flags <- c("No data", "no data")
wet_flags <- c('water','wet','Water','Wet')


assert_that(!sum(ind$L1_tree2 %in% nodata_flags, na.rm = TRUE) &&
            !sum(il$L1_tree2 %in% nodata_flags, na.rm = TRUE) &&
            !sum(det$L1_tree2 %in% nodata_flags, na.rm = TRUE),
            msg = "'No data' found for second tree in IN or IL or Detroit; this case not handled by the code.")


cat("Found ", sum(ind$L1_tree1 %in% wet_flags), " wet corners in Indiana.\n",
    sep = '')
cat("Found ", sum(il$L1_tree1 %in% wet_flags), " wet corners in Illinois.\n",
    sep = '')
cat("Found ", sum(det$L1_tree1 %in% wet_flags), " wet corners in Detroit.\n",
    sep = '')
cat("Found ", sum(ind$L1_tree1 %in% nodata_flags), " 'no data' corners in Indiana.\n",
    sep = '')
cat("Found ", sum(il$L1_tree1 %in% nodata_flags), " 'no data' corners in Illinois.\n",
    sep = '')
cat("Found ", sum(det$L1_tree1 %in% nodata_flags), " 'no data' corners in Detroit.\n",
    sep = '')

## Our density/biomass calculations are on a per-land-area basis, excluding water area
ind <- ind %>% filter(!(L1_tree1 %in% c(nodata_flags, wet_flags)))
il <- il %>% filter(!(L1_tree1 %in% c(nodata_flags, wet_flags)))
det <- det %>% filter(!(L1_tree1 %in% c(nodata_flags, wet_flags)))

cat("Using ", nrow(ind), " corners in Indiana.\n", sep = '')
cat("Using ", nrow(il), " corners in Illinois.\n", sep = '')
cat("Using ", nrow(det), " corners in Detroit.\n", sep = '')

## Converting L1 (survey abbreviation) to L3 (Paleon nomenclature) taxa; currently this overwrites existing L3 in Indiana and Detroit (for Illinois it has already been removed) but ensuring use of current taxon conversion file; in future L3 will not be in the input files and will solely be created here.

warning("waiting on Jody for issue #71 to have 'No tree' in conversion file so NAs not introduced")

## Indiana conversion
spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 2000) %>% 
    filter(domain == indiana_conversion_domain) %>%
    select(level1, level3a)

ind <- ind %>% select(-L3_tree1, -L3_tree2, -L3_tree3, -L3_tree4) %>%  ## so we can replace these columns
    left_join(spec_codes, by = c('L1_tree1' = 'level1')) %>% rename(L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree2' = 'level1')) %>% rename(L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree3' = 'level1')) %>% rename(L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree4' = 'level1')) %>% rename(L3_tree4 = level3a)

assert_that(sum(is.na(ind$L3_tree1)) == sum(is.na(ind$L1_tree1)) &&
   sum(is.na(ind$L3_tree2)) == sum(is.na(ind$L1_tree2)) &&
   sum(is.na(ind$L3_tree3)) == sum(is.na(ind$L1_tree3)) &&
   sum(is.na(ind$L3_tree4)) == sum(is.na(ind$L1_tree4)),
   msg = "Apparently some L1 taxa are missing from the Indiana L1 to L3 conversion table.")

## Illinois conversion
spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 2000) %>% 
    filter(domain == illinois_conversion_domain) %>%
    select(level1, level3a) 

il <- il %>% select(-L3_tree1, -L3_tree2, -L3_tree3, -L3_tree4) %>% 
    left_join(spec_codes, by = c('L1_tree1' = 'level1')) %>% rename(L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree2' = 'level1')) %>% rename(L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree3' = 'level1')) %>% rename(L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree4' = 'level1')) %>% rename(L3_tree4 = level3a)

assert_that(sum(is.na(il$L3_tree1)) == sum(is.na(il$L1_tree1)) &&
   sum(is.na(il$L3_tree2)) == sum(is.na(il$L1_tree2)) &&
   sum(is.na(il$L3_tree3)) == sum(is.na(il$L1_tree3)) &&
   sum(is.na(il$L3_tree4)) == sum(is.na(il$L1_tree4)),
   msg = "Apparently some L1 taxa are missing from the Illinois L1 to L3 conversion table.")

## Detroit conversion
spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 2000) %>% 
    filter(domain == detroit_conversion_domain) %>%
    select(level1, level3a)


det <- det %>% select(-L3_tree1, -L3_tree2, -L3_tree3, -L3_tree4) %>%  ## so we can replace these columns
    left_join(spec_codes, by = c('L1_tree1' = 'level1')) %>% rename(L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree2' = 'level1')) %>% rename(L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree3' = 'level1')) %>% rename(L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('L1_tree4' = 'level1')) %>% rename(L3_tree4 = level3a)

assert_that(sum(is.na(det$L3_tree1)) == sum(is.na(det$L1_tree1)) &&
   sum(is.na(det$L3_tree2)) == sum(is.na(det$L1_tree2)) &&
   sum(is.na(det$L3_tree3)) == sum(is.na(det$L1_tree3)) &&
   sum(is.na(det$L3_tree4)) == sum(is.na(det$L1_tree4)),
   msg = "Apparently some L1 taxa are missing from the Detroit L1 to L3 conversion table.")

## Note that check for NAs appearing in L3 that did not occur in L1 replaces a former
## replacement of NAs with 'No tree', which is not accurate since
## a non-available conversion is not the same as no tree appearing.
## If we really don't know the L3 taxon, we should convert to 'Unknown tree'.

## make sure township names have the state in front of them for later merging of state data:
ind <- ind %>% mutate(twp = paste0('IN_', TRP)) %>% dplyr::select(-TRP)
il <- il %>% mutate(twp = paste0('IL_', TRP)) %>% dplyr::select(-TRP)
det <- det %>% mutate(twp = paste0('MI_', TRP)) %>% dplyr::select(-TRP)

assert_that(all(ind$state == 'IN') && all(il$state == 'IL'),
    msg = "State field missing from one or more rows in IN or IL.")

ind <- ind %>% rename(dist1 = chainstree, dist2 = chainstree2, dist3 = chainstree3, dist4 = chainstree4,
                      diameter1 = diameter, bearing1 = bearing, bearingdir1 = bearingdir, degrees1 = degrees)
il <- il %>% rename(dist1 = chainstree, dist2 = chainstree2, dist3 = chainstree3, dist4 = chainstree4,
                    diameter1 = diameter, bearing1 = bearing, bearingdir1 = bearingdir, degrees1 = degrees)
det <- det %>% rename(dist1 = chainstree, dist2 = chainstree2, dist3 = chainstree3, dist4 = chainstree4,
                    diameter1 = diameter, bearing1 = bearing, bearingdir1 = bearingdir, degrees1 = degrees)

ind$bearing1 <- paste(ind$bearing1, ind$bearingdir1, sep = '_')
ind$bearing2 <- paste(ind$bearing2, ind$bearingdir2, sep = '_')
ind$bearing3 <- paste(ind$bearing3, ind$bearingdir3, sep = '_')
ind$bearing4 <- paste(ind$bearing4, ind$bearingdir4, sep = '_')
il$bearing1 <- paste(il$bearing1, il$bearingdir1, sep = '_')
il$bearing2 <- paste(il$bearing2, il$bearingdir2, sep = '_')
il$bearing3 <- paste(il$bearing3, il$bearingdir3, sep = '_')
il$bearing4 <- paste(il$bearing4, il$bearingdir4, sep = '_')
det$bearing1 <- paste(det$bearing1, det$bearingdir1, sep = '_')
det$bearing2 <- paste(det$bearing2, det$bearingdir2, sep = '_')
det$bearing3 <- paste(det$bearing3, det$bearingdir3, sep = '_')
det$bearing4 <- paste(det$bearing4, det$bearingdir4, sep = '_')

## We have some corners in IN & IL that are missing years. Based on exploratory analysis of their locations
## it is safe to assume these points were surveyed at a similar time as the points around them.

distances <- rdist(ind[ind$year == 9999, c('x', 'y')],
              ind[ind$year != 9999, c('x', 'y')])
closest <- apply(distances, 1, which.min)
ind$year[ind$year == 9999] <- ind$year[ind$year != 9999][closest]
assert_that(!sum(is.na(ind$year)) && min(ind$year) >= 1799 && max(ind$year) <= 1849,
    msg = "Unexpected missing year or year outside 1799-1849 range in IN.")
rm(distances)

distances <- rdist(il[il$year == 9999, c('x', 'y')],
              il[il$year != 9999, c('x', 'y')])
closest <- apply(distances, 1, which.min)
il$year[il$year == 9999] <- il$year[il$year != 9999][closest]
assert_that(!sum(is.na(il$year)) && min(il$year) >= 1799 && max(il$year) <= 1849,
    msg = "Unexpected missing year or year outside 1799-1849 range in IL.")
rm(distances)

## create a survey year variable that corresponds to survey year correction factors
## for IL/IN, '1825+' vs. '< 1825' and state uniquely defines correction factor - don't need region info 
## inil <- inil %>% mutate(surveyyear = ifelse(year >= 1825, '>=1825', '<=1824'))

ind <- ind %>% mutate(surveyyear = ifelse(year >= 1825, '>=1825', '<=1824'))
det <- det %>% mutate(surveyyear = ifelse(year >= 1825, '>=1825', '<=1824'))
assert_that(all(det$surveyyear == "<=1824"), msg = "Found Detroit data with year >1824.")

## As of April 2019, we now split IL into three sets
il <- il %>% mutate(surveyyear = ifelse(year <= 1810, '<=1810',
                                            ifelse(year >= 1838, ">=1838", ">1810  <1838")))


columns_to_keep <- c("x","y","twp","surveyyear",
                     "L1_tree1", "L1_tree2", "L1_tree3", "L1_tree4",
                     "L3_tree1", "L3_tree2", "L3_tree3", "L3_tree4",
                     "bearing1", "bearing2", "bearing3", "bearing4",
                     "degrees1", "degrees2", "degrees3","degrees4",
                     "dist1", "dist2", "dist3", "dist4",
                     "diameter1", "diameter2", "diameter3", "diameter4",
                     "cornerid", "typecorner","state", "point_id", "vegtype")

ind <- ind[columns_to_keep] 
il <- il[columns_to_keep]
det <- det[columns_to_keep]

inildet <- rbind(ind, il, det)

inildet <- inildet %>% rename(diam1 = diameter1, diam2 = diameter2, diam3 = diameter3, diam4 = diameter4)

## Change 88888/99999 to NA but not in in bearing columns because an NA plus a degrees of 0 means a cardinal direction while 88888/99999 is unknown
## 88888/99999 in taxa have been dealt with in L1->L3 conversion
cols <- c("degrees1", "degrees2", "degrees3","degrees4",
          "dist1", "dist2", "dist3", "dist4",
          "diam1", "diam2", "diam3", "diam4")

inildet[ , cols] <- sapply(inildet[ , cols], convert_to_NA, missingCodes = c(88888,99999))

notree <- inildet %>% filter(L1_tree1 %in% c('No tree', 'no tree'))
assert_that(sum(is.na(notree$dist1) & is.na(notree$dist2) & is.na(notree$dist3) & is.na(notree$dist4) &
       is.na(notree$diam1) & is.na(notree$diam2) & is.na(notree$diam3) & is.na(notree$diam4))
       == nrow(notree),
    msg = "Found non-NA distances or diameters for no tree points in IN or IL or Detroit.")


inildet[ , paste0('az', 1:4)] <- get_angle_inil(as.matrix(inildet[ , paste0('bearing', 1:4)]),
                           as.matrix(inildet[ , paste0('degrees', 1:4)]))

assert_that(max(inildet[ , paste0('az', 1:4)], na.rm = TRUE) < 360 &&
            min(inildet[ , paste0('az', 1:4)], na.rm = TRUE) >= 0,
    msg = "Found azimuths outside of 0-359 in IN/IL/Detroit.")


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

inildet <- inildet %>% mutate(sectioncorner = ifelse(inildet$cornerid %in% intsec | inildet$cornerid %in% extsec, 
                                               'section', 'quartersection'), 
                        corner = ifelse(inildet$cornerid %in% intsec | inildet$cornerid %in% intqtr, 
                                        'internal', 'external'))

inildet <- inildet[final_columns]

## ----------------------------------DATA CLEANING: SOUTHERN MI --------------------------------------------------

somi <- read_csv(file.path(raw_data_dir, southern_michigan_file), guess_max = 100000)
somi <- somi %>% mutate(point_id = seq_len(nrow(somi)), vegtype = NA, state = "SoMI")

## Per GH issue #59 we are not using Isle Royale at this point.

if(FALSE) {
    ## Isle Royale data provided separately but in same format as southern MI data
    nomi_extra <- read_csv(file.path(raw_data_dir, northern_michigan_supp_file), guess_max = 100000)
    nomi_extra <- nomi_extra %>% mutate(point_id = seq_len(nrow(nomi_extra)), vegtype = NA, state = "NoMI_extra")
}

## These codes are not used in southern Michigan so don't need to do this filtering:
## somi <- somi %>% filter(!(species1 %in% c('No data', 'Water', 'Wet')))

## formerly we check for species1 and species2 being missing but having a positive diameter but
## there is only one case of missing species and existing diameter and that is 4th tree in
## case where it is not among the closest two. 

##  converting level 1 species to level 3 species

## actual so MI points
spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 2000) %>% 
    filter(domain == southern_michigan_conversion_domain) %>%
    select(level1, level3a) 

## just call this 'somi' if/when take out code for nomi_extra
somi_lower <- somi %>% 
    left_join(spec_codes, by = c('species1' = 'level1')) %>% rename(L1_tree1 = species1, L3_tree1 = level3a) %>%
    left_join(spec_codes, by = c('species2' = 'level1')) %>% rename(L1_tree2 = species2, L3_tree2 = level3a) %>%
    left_join(spec_codes, by = c('species3' = 'level1')) %>% rename(L1_tree3 = species3, L3_tree3 = level3a) %>%
    left_join(spec_codes, by = c('species4' = 'level1')) %>% rename(L1_tree4 = species4, L3_tree4 = level3a)

if(FALSE) {
    ## Isle Royale points
    spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 1000) %>% 
        filter(domain == upper_midwest_conversion_domain) %>%
        select(level1, level3a) 
    
    nomi_extra <- nomi_extra %>% 
        left_join(spec_codes, by = c('species1' = 'level1')) %>% rename(L1_tree1 = species1, L3_tree1 = level3a) %>%
        left_join(spec_codes, by = c('species2' = 'level1')) %>% rename(L1_tree2 = species2, L3_tree2 = level3a) %>%
        left_join(spec_codes, by = c('species3' = 'level1')) %>% rename(L1_tree3 = species3, L3_tree3 = level3a) %>%
        left_join(spec_codes, by = c('species4' = 'level1')) %>% rename(L1_tree4 = species4, L3_tree4 = level3a)

    ## Remove the entries with sec_corner == "Check  -  Double" on Isle Royale
    ## as these have two points at a location. Keep the entries with 'Check'
    ## as Charlie Cogbill checked these; they are valid data but location is a
    ## bit off and puts them in the Lake.
    nomi_extra <- nomi_extra %>% filter(!sec_corner %in% c("Check - double"))

    ## of course 'somi' is a misnomer given presence of northern MI data
    somi <- rbind(somi_lower, nomi_extra)
} else {
    somi <- somi_lower
}

if(sum(is.na(somi$L3_tree1)) != sum(is.na(somi$L1_tree1)) ||
   sum(is.na(somi$L3_tree2)) != sum(is.na(somi$L1_tree2)) ||
   sum(is.na(somi$L3_tree3)) != sum(is.na(somi$L1_tree3)) ||
   sum(is.na(somi$L3_tree4)) != sum(is.na(somi$L1_tree4)))
    cat("Apparently some L1 taxa are missing from the southern Michigan L1 to L3 conversion table.\n")

## determine subdomain for use with correction factors
assert_that(nrow(somi) == length(grep('[EW]', somi$range)),
    msg = "Can't assign surveyyear for some southern Michigan sites.")
## shorthand (and unique) for 'SE - E of central Meridian S of tension'
surveyyear <- rep('<=1824', nrow(somi))  
## shorthand (and unique) for 'SW - W of central Meridian S of tension'
surveyyear[grep('W', somi$range)] <- '1825 -1835'
## these townships in northern part of southern Michigan have their own correction factors
special <- read_csv(file.path(conversions_data_dir, michigan_special_township_file))
surveyyear[somi$twnrng %in% special$TWNRNG] <- ">=1831"
if(FALSE) {
    ## Isle Royale in UP: 1847-1848 is shorthand (and unique) for 'Isle Royale, 1847-1848'
    surveyyear[somi$point_y > 900000] <- '1847-1848'
}

somi$surveyyear <- surveyyear

## make sure township names have the state in front of them:
somi <- somi %>% mutate(twp = paste0('MI_', twnrng)) 

## The az_360 columns were calculated using the quadrant and the az_x (x = 1:4), values.
## e.g., if the quadrant number was 1, the AZ as read from the mylar maps was used as is.
## If the quadrant number was 2, az1_360 is 180-az1.
## Quadrants of 3 were 180+az1 and quadrants of 4 were 360-az1.
somi <- somi %>% select(-az1, -az2, -az3, -az4) %>% 
    rename(az1 = az1_360, az2 = az2_360, az3 = az3_360, az4 = az4_360)  

somi <- somi %>% mutate(corner = ifelse(sec_corner == "Extsec", 'external', 'internal'),
                    sectioncorner = ifelse(cornertype == 'section',  'section', 'quartersection'))

somi <- somi %>% rename(x = point_x, y = point_y)

## Per issue 49, quarter-section points in southern Michigan are undersampled and sampling
## appears to vary with forest type/density, so removing them.

count <- somi %>% filter(sectioncorner == 'quartersection') %>% summarize(count = n())
cat("Removing all", unlist(count), "quarter section points from southern Michigan outside Detroit.\n")
somi <- somi %>% filter(sectioncorner == 'section')

cat("Using ", nrow(somi), " corners in southern Michigan (not counting Detroit).\n", sep = '')

somi <- somi[final_columns]

## ----------------------------------DATA CLEANING: Upper Midwest -------------------------------------------------

## Data cleaning modified from Simon Goring's witness tree code
## (https://github.com/PalEON-Project/WitnessTrees/blob/master/R/process_raw/step.one.clean.bind_v1.4.R)

##  Original warning from Simon Goring's processing:
## 
##  Binding and cleaning the Wisconsin and Minnesota Public Lands Surveys, data is
##  sourced from the Mladenoff Lab at the University of Wisconsin.  The lab has
##  granted us permission to use the data, but not permission to distribute the
##  original datasets.  A version of the Minnesota data can be obtained from 
##  http://deli.dnr.state.mn.us/metadata/pveg_btreept3.html
##  The wisconsin data may be obtained by contacting David Mladenoff at:
##  mladenoff@wisc.edu

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
mn <- mn %>% mutate(point_id = seq_len(nrow(mn)))
nomi <- nomi %>% mutate(point_id = seq_len(nrow(nomi)), vegtype = NA)

names(wi) <- tolower(names(wi))
names(mn) <- tolower(names(mn))
names(nomi) <- tolower(names(nomi))

wi <- wi %>% rename(year = year_)

nomi <- nomi %>% rename(diam1 = dbh1, diam2 = dbh2, diam3 = dbh3, diam4 = dbh4,
                        sp1 = spp1, sp2 = spp2, sp3 = spp3, sp4 = spp4,
                        az1 = azimuth, az2 = azimuth2, az3 = azimuth3, az4 = azimuth4)

assert_that(all(wi$rangdir %in% c(0, 2, 4)), 
    msg = "Unexpected 'rangdir' values found in Wisconsin.")

##  The wisconsin Range is set as a single value, the 'E' and 'W' codes are in
##  RANGDIR.  Looking at the data it also looks like there are a few ranges
##  that are miscoded with slight differences from cases Simon was checking. 
wi <- wi %>% filter(rangdir != 0) %>%   ## Single point that seems to have no data associated with it, so remove it.
    mutate(rangdir = ifelse(rangdir == '2', 'W', 'E')) %>%
    mutate(rangdir = ifelse(rangdir == 'W' & y_alb > 1050000 & x_alb > 4.6e5, 'E', rangdir)) %>%
    mutate(rangdir = ifelse(rangdir == 'E' & y_alb > 1050000 & x_alb < 4.3e5, 'W', rangdir)) %>%
    mutate(rng = paste0(as.character(range), rangdir)) %>%
    select(-range, -rangdir)

## per GH issue #60, 'point' field is now available as 'pnt'
nomi <- nomi %>% rename(point = pnt)

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

assert_that(sum(mn$vegtype %in% waterTypes & mn$sp1 == '_' & (mn$sp2 != '_' | mn$sp3 != '_' | mn$sp4 != '_')) == 0,
            msg = "MN water corners have trees with species info.")

## '_' values in areas with water tagged as QQ to follow WI labelling
mn <- mn %>% mutate(sp1 = ifelse(vegtype %in% waterTypes & sp1 == '_', 'QQ', sp1),
                    sp2 = ifelse(vegtype %in% waterTypes & sp2 == '_', 'QQ', sp2),
                    sp3 = ifelse(vegtype %in% waterTypes & sp3 == '_', 'QQ', sp3),
                    sp4 = ifelse(vegtype %in% waterTypes & sp4 == '_', 'QQ', sp4))

## other '_' values treated as NA; a very small number of these having non-missing distance or diameter
mn <- mn %>% mutate(sp1 = convert_to_NA(sp1, '_'),
              sp2 = convert_to_NA(sp2, '_'),
              sp3 = convert_to_NA(sp3, '_'),
              sp4 = convert_to_NA(sp4, '_'))

## exclude water points -- all water as well as points with 1 tree are remaining water, per issue #35
numQQ <- apply(mn[ , paste0('sp', 1:4)], 1, function(x) sum(x == 'QQ', na.rm = TRUE))
cat("Found ", sum(numQQ > 2), " wet corners in Minnesota.\n", sep = '')
mn <- mn %>% filter(numQQ <= 2)

cat("Keeping ", sum(numQQ %in% c(1,2)), " corners in possibly wet areas with at least two trees.\n")

## Next exclude no-tree points with unclear vegtype values

## these points occur mostly in northern Minnesota (Boundary Waters) or
## are in east-west straight lines, suggesting not usable
nr <- nrow(mn)
mn <- mn %>% filter(!(is.na(sp1) & is.na(sp2) & is.na(sp3) & is.na(sp4) & vegtype == '_'))
cat("Excluding ", nr - nrow(mn), " points with missing taxa and missing ecotype in Minnesota.\n")

## forest, grove, bottom, pine grove seem inconsistent with lack of trees, so exclude these points
forestedTypes <- c('F', 'G', 'H', 'J')
nr <- nrow(mn)
mn <- mn %>% filter(!(is.na(sp1) & is.na(sp2) & is.na(sp3) & is.na(sp4) & vegtype %in% forestedTypes))
cat("Excluding ", nr - nrow(mn), " points with missing taxon in forested areas in Minnesota.\n")

mn_survey <- read_csv(file.path(raw_data_dir, minnesota_survey_file), guess_max = 10000)
mn_survey <- mn_survey %>% mutate(TOWN = paste0('T', formatC(TOWN, width=3, flag='0'), 'N'),
                                  RANG = paste0('R', formatC(RANG, width=2, flag='0'), RDIR)) %>%
    dplyr::select(TOWN, RANG, YEAR)
names(mn_survey) <- tolower(names(mn_survey))

mn <- mn %>% left_join(mn_survey, by = c('twp' = 'town', 'rng' = 'rang'))

## two points with twp/rng not in survey file so need to estimate year based on closest points
distances <- rdist(mn[is.na(mn$year), c('x_alb','y_alb')],
                   mn[!is.na(mn$year), c('x_alb','y_alb')])
closest <- apply(distances, 1, which.min)
mn$year[is.na(mn$year)] <- mn$year[!is.na(mn$year)][closest]

assert_that(sum(is.na(mn$year)) == 0 && min(mn$year) >= 1847 && max(mn$year) <= 1907,
    msg = "Unexpected missing year or year outside 1847-1907 range.")

## '<=1853' vs. '>=1854' plus state uniquely defines correction factors without need for region info
mn <- mn %>% mutate(surveyyear = ifelse(year <= 1853, "<=1853", ">=1854"))

## requires about 1 GB RAM
cat("Found ", sum(wi$year == 0), " points without year in Wisconsin.\n")
distances <- rdist(wi[wi$year == 0, c('x_alb','y_alb')],
                   wi[wi$year != 0, c('x_alb','y_alb')])
closest <- apply(distances, 1, which.min)
wi$year[wi$year == 0] <- wi$year[wi$year != 0][closest]
rm(distances)

assert_that(sum(is.na(wi$year)) == 0 && min(wi$year) >= 1832 && max(wi$year) <= 1891,
    msg = "Unexpected missing year or year outside 1847-1907 range.")

## '<=1845' vs. '>=1846' plus state uniquely defines correction factors without need for region info
wi <- wi %>% mutate(surveyyear = ifelse(year <= 1845, "<=1845", ">=1846"))

## In essentially all of these cases, no info on trees 2-4,
## (and in remaining cases, have '0', 'NL' or 'NO' indicating no trees/information)
## so assume these are fully water
notree_vals <- c('0','NO','NL')
assert_that(sum(is.na(wi$sp2[wi$sp1 == 'QQ'])) + sum(wi$sp2[wi$sp1 == 'QQ'] %in% notree_vals) ==
            sum(wi$sp1 == 'QQ'), msg = "Found tree taxa in Wisconsin water points.")
assert_that(sum(is.na(wi$sp3[wi$sp1 == 'QQ'])) + sum(wi$sp3[wi$sp1 == 'QQ'] %in% notree_vals) ==
            sum(wi$sp1 == 'QQ'), msg = "Found tree taxa in Wisconsin water points.")
assert_that(sum(is.na(wi$sp4[wi$sp1 == 'QQ'])) + sum(wi$sp4[wi$sp1 == 'QQ'] %in% notree_vals) ==
            sum(wi$sp1 == 'QQ'), msg = "Found tree taxa in Wisconsin water points.")
cat("Found ", sum(wi$sp1 == 'QQ'), " water points in Wisconsin.\n")
wi <- wi %>% filter(sp1 != 'QQ')

## based on parsing notes field, we do include some points that seem to be no-tree points

## (non-water and non-"tree is point/corner/post")
miss <- nomi %>% filter(is.na(sp1) & is.na(sp2) & is.na(sp3) & is.na(sp4))
nomi <- nomi %>% filter(!(is.na(sp1) & is.na(sp2) & is.na(sp3) & is.na(sp4)))

cat("Attempting to determine status of ", nrow(miss), " points in northern Michigan with all taxa missing.\n")

## points without any taxa and no notes are ambiguous so remove
cat("Omitting ", sum(is.na(miss$notes)), " points in northern Michigan with no notes.\n")
miss <- miss %>% filter(!is.na(notes))
## water points
waterText <- c("(MARSH|POND|LAKE|WATER|RIVER|SWAMP|BROOK|STREAM|LK MICHIGAN)")
cat("Omitting ", sum(grepl(waterText, miss$notes, ignore.case = TRUE)), " points in northern Michigan with indications of water and no taxa info.\n")
miss <- miss %>% filter(!grepl(waterText, notes, ignore.case = TRUE))
## points where tree is the corner so ambiguous what density would be
## some of these might be trees at previous points, but can't determine which
treeAsPostText <- c("(IS POST|IS CORNER|IS QUARTER|IS SECTION|AS CORNER|UPON CORNER|AS 1/4|IS 1/4|FOR 1/4|FOR CORNER|IS PSOT|ISPOST|AS POST|IS A QUARTER|IN CORNER|IS THE QUARTER CORNER|IS A QAURTER CORNER|IS A SECTION|A QUARTER CORNER|IS THE SECTION CORNER|FOR QUARTER|IS THE QAURTER CORNER)")
cat("Omitting ", sum(grepl(treeAsPostText, miss$notes, ignore.case = TRUE)), " points in northern Michigan with corner as post and no taxa info.\n")
miss <- miss %>% filter(!grepl(treeAsPostText, notes, ignore.case = TRUE))
## indications that data lost or not noted
missingDataText <- c("(THEN LOST|INFORMATION NOT GIVEN|NOT IN NOTES|NO INFORMATION|RANDOM NOTES|OMITTED|OMMITTED|CUT OFF|MICROFILM|TREE IS DEAD|TA, RP, PS, WP|PS, RP|TRAIL COURSE)")
cat("Omitting ", sum(grepl(missingDataText, miss$notes, ignore.case = TRUE)), " points in northern Michigan with corner as post and no taxa info.\n")
miss <- miss %>% filter(!grepl(missingDataText, notes, ignore.case = TRUE))
## Charlie's checking indicates these are mostly cases where tree was the corner but surveyor didn't bother to note other trees; density here is ambiguous
noOtherTreeText <- c("NO OTHER")
cat("Omitting ", sum(grepl(noOtherTreeText, miss$notes, ignore.case = TRUE)), " points in northern Michigan marked 'NO OTHER'.\n")
miss <- miss %>% filter(!grepl(noOtherTreeText, notes, ignore.case = TRUE))

## otherwise, notes generally say 'no witness trees', 'no trees convenient', 'no bearing trees', 'no other tree data', 'no trees'; assumed to indicate no-tree points
nomi <- rbind(nomi, miss)

## >=1840 is shorthand (and unique) for 'UP, >=1840'
## >=1836 is shorthand (and unique) for 'north Lower, >=1836'
up_lp <- read_csv(file.path(raw_data_dir, michigan_up_lp_file)) %>% select(FID_, Location)
nomi <- nomi %>% left_join(up_lp, by = c('fid' = 'FID_'))
## Per Charlie Cogbill email, some UP points were surveyed before 1840 but this is really
## just an ad hoc approach to get the right correction factors.
nomi <- nomi %>% mutate(surveyyear = ifelse(Location == "UP", '>=1840', '>=1836'))

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

cat("Using ", nrow(mn), " corners in Minnesota.\n", sep = '')
cat("Using ", nrow(wi), " corners in Wisconsin, but some will excluded later as water points.\n", sep = '')
cat("Using ", nrow(nomi), " corners in northern Michigan.\n", sep = '')


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

spec_codes <- read_csv(file.path(conversions_data_dir, taxa_conversion_file), guess_max = 2000) %>% 
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

## These are points in WI marked as XA/XB/XC where no corner data present.
miss <- umw %>% filter(umw$L3_tree1 %in% nodata_flags)
nr <- nrow(miss)
assert_that(sum(is.na(miss$L3_tree2) | miss$L3_tree2 == "No tree") == nr &&
            sum(is.na(miss$L3_tree3) | miss$L3_tree3 == "No tree") == nr &&
            sum(is.na(miss$L3_tree4) | miss$L3_tree4 == "No tree") == nr, msg = "Strange missing corners in Wisconsin")
cat("Excluding ", nr, " points with no data in Wisconsin.\n")
umw <- umw %>% filter(!umw$L3_tree1 %in% nodata_flags)

#  Points within a township are either sections or quartersections.  This
#  is the list of points that are sections.  All others are quarter-sections.
sections <- c(2, 5, 8, 11, 14, 18, 21, 24, 27, 30,
              34, 37, 40, 43, 46, 50, 53, 56, 59, 62,
              66, 70, 74, 78, 82,
              87, 89, 91, 93, 95, 98, 100, 102, 104, 106, 108,
              109, 111, 113, 115, 117, 119, 122, 123, 124, 125, 126)

#  These are the points on the outside of each township.
external <- c(109:120, 97:108, 86:96, 122:126)

## illegitimate point values, preventing determination of correction factors
assert_that(all(umw$point %in% 1:126),
    msg = "upper Midwest point values outside 1:126 found.")
umw <- umw %>% filter(point %in% 1:126)

umw <- umw %>% mutate(corner = ifelse(point %in% external, 'external', 'internal'),
                      sectioncorner = ifelse(point %in% sections, 'section', 'quartersection'))

umw <- umw[final_columns]

## --------------- Combine regions and do further cleaning ------------------------------

mw <- rbind(umw, inildet, somi)

##  At this point we need to make sure that the species are ordered by distance
##  so that trees one and two are actually the closest two trees.

taxa <- c(mw$L3_tree1, mw$L3_tree2, mw$L3_tree3, mw$L3_tree4)
assert_that(sum(taxa %in% nodata_flags) == 0, msg = "'No data' points still in combined data.")

## 'Missing' means either dead or the XX flag, which is indeterminate.
## Almost all points with dead trees do not have data on two live trees.
## There are relatively few points with dead trees, so just exclude them.
## It's hard to use filter() in this case becauase want to keep all the NAs.
nr <- nrow(mw)

wh <- mw$L3_tree1 == "Missing" | mw$L3_tree2 == "Missing" |
    mw$L3_tree3 == "Missing" | mw$L3_tree4 == "Missing"
wh[is.na(wh)] <- FALSE
mw <- mw[!wh, ]
cat("Removing ", nr-nrow(mw), " points with any 'Dead' or indeterminable trees as these points generally don't have two live trees for calculation.\n")

nontree_codes <- c("Water", "No tree")

## set dists and azimuths to NA when there is not a tree there (NA, water, no tree) so
## zeroes are not interpreted as 0 distance or as being in same quad

mw <- mw %>% mutate(dist1 = ifelse(is.na(L3_tree1) | L3_tree1 %in% nontree_codes, NA, dist1),
                    dist2 = ifelse(is.na(L3_tree2) | L3_tree2 %in% nontree_codes, NA, dist2),
                    dist3 = ifelse(is.na(L3_tree3) | L3_tree3 %in% nontree_codes, NA, dist3),
                    dist4 = ifelse(is.na(L3_tree4) | L3_tree4 %in% nontree_codes, NA, dist4))

mw <- mw %>% mutate(az1 = ifelse(is.na(L3_tree1) | L3_tree1 %in% nontree_codes, NA, az1),
                    az2 = ifelse(is.na(L3_tree2) | L3_tree2 %in% nontree_codes, NA, az2),
                    az3 = ifelse(is.na(L3_tree3) | L3_tree3 %in% nontree_codes, NA, az3),
                    az4 = ifelse(is.na(L3_tree4) | L3_tree4 %in% nontree_codes, NA, az4))

## per issue #39 we are calculating density and biomass based on trees below veil line
## with correction to get density for trees above the veil line
if(FALSE){
    ## set small trees (below 8 inch veil line) dists to Inf so we find the bigger trees as closest two
    mw <- mw %>% mutate(diam1 = ifelse(diam1 < diameter_cutoff_inches, NA, diam1), 
                        diam2 = ifelse(diam2 < diameter_cutoff_inches, NA, diam2),
                        diam3 = ifelse(diam3 < diameter_cutoff_inches, NA, diam3),
                        diam4 = ifelse(diam4 < diameter_cutoff_inches, NA, diam4))
    
    mw <- mw %>% mutate(dist1 = ifelse(diam1 < diameter_cutoff_inches, Inf, dist1), 
                        dist2 = ifelse(diam2 < diameter_cutoff_inches, Inf, dist2),
                        dist3 = ifelse(diam3 < diameter_cutoff_inches, Inf, dist3),
                        dist4 = ifelse(diam4 < diameter_cutoff_inches, Inf, dist4))
}

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

assert_that(all(mw$dist1 <= mw$dist2, na.rm = TRUE) &&
            all(mw$dist2 <= mw$dist3, na.rm = TRUE) &&
            all(mw$dist3 <= mw$dist4, na.rm = TRUE),
            msg = "Distances not reordered correctly.")

## determine Pair (points where only two trees surveyed) vs 2nQ (four trees surveyed)
## used for correction factors (see issue #42)
## it would be rare to only find two trees if looking for four,
## and if we have zero or one tree, we don't use correction factors anyway
tmp <- as.matrix(mw[ , paste0('L3_tree', 1:4)])
assert_that(length(unique(c(tmp))) == 38,
            msg = "Unexpected level 3a taxa found")

tmp[tmp %in% c('No tree', 'Water')] <- NA
ntree <- apply(tmp, 1,function(x) sum(!is.na(x)))
mw <- mw %>% mutate(point = ifelse(ntree > 2, '2nQ', 'Pair'))


## before additional cleaning of northern Michigan datamost there were decimal distances
## that might be chains. After cleaning, only ~230 decimal distances and all but one are >= 1.5
## and are in the form x.5, so assuming they are links.
if(FALSE) {
    mw <- mw %>% mutate(dist1 = convert_chains_to_links(dist1),
                        dist2 = convert_chains_to_links(dist2),
                        dist3 = convert_chains_to_links(dist3),
                        dist4 = convert_chains_to_links(dist4))
}

save(mw, file = file.path(interim_results_dir, 'cleaned_point.Rda'))

