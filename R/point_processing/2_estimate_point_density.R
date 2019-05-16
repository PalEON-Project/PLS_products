## Estimate tree density at all points, both density of trees above diameter cutoff
## (for density product) and total density (for use in biomass calculation that
## includes trees below the cutoff).

## Run time: approximately 15 seconds

library(readr)
library(dplyr)
library(assertthat)

load(file.path(interim_results_dir, 'cleaned_point.Rda'))

## determine number of trees per point

wet_flags <- c("wet", "Wet")
assert_that(sum(mw$L3_tree1 %in% wet_flags, na.rm = TRUE) + 
            sum(mw$L3_tree2 %in% wet_flags, na.rm = TRUE) == 0, msg = "'wet' found as tree species.")

num_trees <- rep(2, nrow(mw))
num_trees[is.na(mw$L3_tree1) | mw$L3_tree1 == "No tree"] <- 0
num_trees[!is.na(mw$L3_tree1) & mw$L3_tree1 != "No tree" & (mw$L3_tree2 == "No tree" | is.na(mw$L3_tree2))] <- 1

mw <- mw %>% mutate(num_trees = num_trees)

## Apply various exclusion criteria

## remove 1-tree 0-distance (or unknown distance) points as unclear what to use for density
## there are clusters of such points in LP of MI, southern IL and southern IN, and around Green Bay
nr <- nrow(mw)
mw <- mw %>% filter(!(num_trees == 1 & (is.na(dist1) | dist1 == 0)))
cat("Excluded ", nr - nrow(mw), " points with one tree at zero or missing distance.\n")

## can't calculate density for points with missing distances
nr <- nrow(mw)
mw <- mw %>% filter(!(num_trees == 2 & (is.na(dist1) | is.na(dist2))))
cat("Excluded ", nr - nrow(mw), " points with two trees and one or more missing distances.\n")

## can't calculate density for points with two zero distances: 
nr <- nrow(mw)
mw <- mw %>% filter(!(num_trees == 2 & dist1 == 0 & dist2 == 0))
cat("Excluded ", nr - nrow(mw), " points with two trees and two distances equal zero.\n")

## Remove one-tree WI points with indications of water;
## doing this here as it's simplest to make use of the num_trees field
## assuming same vegtype codes as in MN, we should exclude L,M,R,S
waterTypes <- c('L', 'M', 'S', 'R', 'A')
nr <- nrow(mw)
mw <- mw %>% filter(!(state == 'WI' & num_trees %in% c(0,1) & vegtype %in% waterTypes)) %>%
    select(-vegtype)
cat("Excluded ", nr - nrow(mw), " Wisconsin points with zero or one tree in water areas.\n")

## Most of these az=360 values are in southern MI and Charlie has confirmed that this is as expected.
if(any(mw[ , paste0('az', 1:4)] == 360)) {
    mw <- mw %>% mutate(az1 = ifelse(az1 == 360, 0, az1),
                        az2 = ifelse(az2 == 360, 0, az2),
                        az3 = ifelse(az3 == 360, 0, az3),
                        az4 = ifelse(az4 == 360, 0, az4))
}

## Charlie Cogbill has checked all MI trees > 50 inches.
## We should have a limited number of trees > 50 remaining.
## Note allometries mostly wouldn't go above 80 cm = 32 inches.
tmp <- mw %>% filter(diam1 >= 50 | diam2 >= 50)
assert_that(nrow(tmp) < 500,
            msg = "more than 500 points with trees with diameter greater than 50 inches")

## keep 2-tree points regardless of distances and truncate density (at say 1000 for now)
## use 2-tree points with small or NA diameter trees for density calculation,

## Treat all 2nQ points in southern Michigan as Pair as that is what we think
## the sampling scheme was and we don't have correction factors for 2nQ in this area.
assert_that(nrow(mw[mw$state == "SoMI" & mw$point == "2nQ",]) < 50,
            msg = "more than 50 2nQ points in southern Michigan")
mw <- mw %>% mutate(point = ifelse(state == "SoMI" & point == "2nQ", "Pair", point))

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

print(table(mw$state[!is.na(mw$density)]))

save(mw, file = file.path(interim_results_dir, 'point_with_density.Rda'))



