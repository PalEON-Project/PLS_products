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

## TODO: fix 2nQ vs P
mw <- mw %>% mutate(density = calc_stem_density(mw, corr_factors))

save(mw, file = 'point_with_density.Rda')

## TODO: look for and perhaps  truncate anomalously high stem density


