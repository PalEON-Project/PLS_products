## Estimate tree density at all points, both density of trees above diameter cutoff
## (for density product) and total density (for use in biomass calculation that
## includes trees below the cutoff).

## Run time: approximately 15 seconds

library(readr)
library(dplyr)

load(file.path(interim_results_dir, 'cleaned_point.Rda'))

## determine number of trees per point

num_trees <- rep(2, nrow(mw))
## only 1 NA in L3_tree1; check back on where it came from
num_trees[is.na(mw$L3_tree1) | mw$L3_tree1 == "No tree"] <- 0
num_trees[!is.na(mw$L3_tree1) & mw$L3_tree1 != "No tree" & (mw$L3_tree2 == "No tree" | is.na(mw$L3_tree2))] <- 1

mw <- mw %>% mutate(num_trees = num_trees)

## Apply various exclusion criteria

## remove 1-tree 0-distance (or unknown distance) points as unclear what to use for density
## there are clusters of such points in LP of MI, southern IL and southern IN, and around Green Bay
## 3413 points
mw <- mw %>% filter(!(num_trees == 1 & (is.na(dist1) | dist1 == 0)))

## can't calculate density for points with missing distances: 1683 points
mw <- mw %>% filter(!(num_trees == 2 & (is.na(dist1) | is.na(dist2))))

## remove points where both distances are less than 1 -- approximately 240 points of which
## ~150 are in four cells in an east-west row in northern lower Michigan;
## presumably a decimal point transciption error
mw <- mw %>% filter(!(!is.na(dist1) & dist1 > 0 & dist1 < 1 &
                      !is.na(dist2) & dist2 > 0 & dist2 < 1))


## Remove one-tree WI points with indications of water;
## doing this here as it's simplest to make use of the num_trees field
## assuming same vegtype codes as in MN, we should exclude L,M,R,S
## however there are codes that are the numbers 2,3,4,5,7,8; not sure what these mean
waterTypes <- c('L', 'M', 'S', 'R', 'A')
mw <- mw %>% filter(!(state == 'WI' & num_trees == 1 & vegtype %in% waterTypes)) %>%
    select(-vegtype)

## Most of these az=360 values are in so MI and given that Charlie cleaned the azimuths there,
## these should be fine.
if(any(mw[ , paste0('az', 1:4)] == 360)) {
    warning("Found some azimuths = 360")
    mw <- mw %>% mutate(az1 = ifelse(az1 == 360, 0, az1),
                        az2 = ifelse(az2 == 360, 0, az2),
                        az3 = ifelse(az3 == 360, 0, az3),
                        az4 = ifelse(az4 == 360, 0, az4))
}


## about 40 trees > 100 in diameter; omit these points as likely erroneous and would induce huge biomass
## note allometries mostly wouldn't go above 80 cm = 32 inches
mw <- mw %>% filter(!( (!is.na(diam1) & diam1 >= max_diam_inch) | (!is.na(diam2) & diam2 > max_diam_inch) ))

## keep 2-tree points regardless of distances and truncate density (at say 1000 for now)
## use 2-tree points with small or NA diameter trees for density calculation,

## keep 1-tree points and set density to something like 1 tree per unit circle of radius (150 links for now)
## ~3.5 trees/ha

# about 10k points with points known to be in same azimuth; about 11k with missing az1, 8500 missing az2

## Sensitivity analyses:
## w/ and w/o 1-tree points - some increase w/o 1-tree points but seems not generally more than 15-20% per cell

## Correction factors to account for surveyor sampling 'design'
corr_factors <- read_csv(file.path(conversions_data_dir, correction_factors_file),
                         na = c("", "NA", "na"))

## fix state names to match tree data
names_df <- data.frame(state = c('Indiana','Illinois','Michigan','Minnesota','Wisconsin'),
                       new_state = c('IN','IL','MI','MN','WI'), stringsAsFactors = FALSE)
corr_factors <- corr_factors %>% left_join(names_df, by = c('state' = 'state')) %>%
    dplyr::select(-state) %>% rename(state = new_state)

## per issue #39, we will omit phi (veil line) correction when computing point biomass
## as we include biomass of trees below the veil line
## therefore compute two density estimates here
mw <- mw %>% mutate(density = calc_stem_density(mw, corr_factors),
                    density_for_biomass = calc_stem_density(mw, corr_factors, use_phi = FALSE))

save(mw, file = file.path(interim_results_dir, 'point_with_density.Rda'))



