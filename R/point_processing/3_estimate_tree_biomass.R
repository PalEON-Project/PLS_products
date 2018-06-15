## Estimate tree-level biomass using diameter of trees and PEcAn allometries

## Run-time: approximately 40 minutes when using allometry parameters shared by cell

library(dplyr)
library(readr)
library(PEcAn.allometry)
library(PEcAn.logger)

cols_conv <- cols(
  level3a = col_character(),
  pecan_allometry_spcd = col_character(),  ## multiple spcd separated by ; for non-exact matches
  pecan_allometry_common_names = col_character()
)


load(file.path(interim_results_dir, 'point_with_density.Rda'))

## MASS masks select, so need dplyr::select (MASS -- comes from Pecan)

taxa_conversion <- read_csv(file.path(conversions_data_dir, pls_to_pecan_conversion_file),
                            col_types = cols_conv) %>%
    dplyr::select(level3a, pecan_allometry_spcd)

## only taxa with a PEcAn allometry (hopefully this will include all trees in dataset)
taxa_conversion <- taxa_conversion %>% 
    mutate(pecan_allometry_spcd = gsub(";", ",", pecan_allometry_spcd))

mw <- mw %>% left_join(taxa_conversion, by = c("L3_tree1" = "level3a")) %>%
    rename(pecan1 = pecan_allometry_spcd) %>%
    left_join(taxa_conversion, by = c("L3_tree2" = "level3a")) %>%
    rename(pecan2 = pecan_allometry_spcd) 

## based on level3a_to_pecan_v0.2.csv
## not fit: 491 (dogwood), 540 (ash), 701 (ironwood)
## one allometry: 970 (elm), 920 (willow), 762 (cherry), 332 (buckeye)
## two allometries: 901 (locust), 731 (sycamore), 241 (cedar/juniper)
## three allometries: 951 (basswood), 351 (seems very heterogeneous) (alder), 71 (perhaps fine) (tamarack)

## 97,94,95 (spruce) had to drop an allometry and wide scatter but a lot of data and looks same as when fit for FIA

## based on FIA:
## 970: 318;802;541;731  (few elm allometries)
## 762: 318;802;541;731  (few small cherries in modern times, so assume black cherry)
## 491, 701: 315;471;491;319;355;356;357;931;935;391;701;660;761;500;501;502;763;662;663;765;766
## 540: 541;543;544
## 332, 901, 920, 241, 951, 351, 71 seemingly used
## 731 seemingly used (as in FIA) but could use 318;802;541;731

## updates made to produce level3a_to_pecan_v0.3.csv

## Treat Unknown tree as generic hardwood.
## Even if we want to omit Unknown tree, we need a placeholder so that NA
## (from Unknown tree) values don't cause problems in allometry fitting
mw <- mw %>% mutate(pecan1 = ifelse(is.na(pecan1), '318,802,541,731', pecan1),
                    pecan2 = ifelse(is.na(pecan2), '318,802,541,731', pecan2))

## Fit allometry models for all taxa
unique_pecan_allom <- unique(c(mw$pecan1, mw$pecan2))
pecan_taxa <- lapply(seq_len(length(unique_pecan_allom)), function(i) 
    data.frame(spcd = as.numeric(strsplit(unique_pecan_allom[i], split = ',')[[1]])))
names(pecan_taxa) <- unique_pecan_allom

## first check which allometries are already in the allom_dir and don't redo
allom_stats = load.allom(allom_dir)

## create via loop rather than passing list of dataframes because of PEcAn bug 
for(i in seq_along(pecan_taxa)) {
    if(!names(pecan_taxa)[i] %in% names(allom_stats)) {
        cat("Fitting ", names(pecan_taxa)[i], ".\n")
        allom_stats[[names(pecan_taxa)[i]]] <- try(AllomAve(pecan_taxa[i], ngibbs = 1000, components = 6,
                                                            outdir = allom_dir, dmin = 10, dmax = 150))
    }}
##allom_stats = AllomAve(pecan_taxa, ngibbs=1000, components = 6, outdir = allom_dir, dmin = 10, dmax = 150)


## site effects seem to product crazy allometries - mu0 and mu1 have outliers and taus can be big; how many allometries do we need for stability?

## We want to run the allometry prediction for all trees of a given taxon at a given point simultaneously
## because this allows for shared allometric parameters (but different tree effects) for trees of the same taxon.
## In this case we simply run the prediction for all trees at a point at once.

## Actually, arguably we should use shared parameters for all trees in a given cell,
## but hard to know what spatial scale to assume homogeneity at.
## Given the increased computational time of running separately for each point,
## we are currently (based on shared_params_in_cell) doing for all trees in a cell.

mw <- mw %>% mutate(biomass1 = NA, biomass2 = NA)

## This takes ?? (3.5) hours on smeagol on one core.
## Parallelizing would require data output manipulation because we write into rows of 'mw'
## Note that every call loads in the allometry material for all of the fitted allometries,
## so it would be easy to make this more efficient.
set.seed(1)

## NAs in diam become in NAs in biomass (with warnings)
if(!shared_params_in_cell) {
    non_empty_points <- seq_len(nrow(mw))
    non_empty_points <- non_empty_points[mw$num_trees > 0]
    for(pnt in non_empty_points) {
        if(mw$num_trees[pnt] == 2) {
            pred <- allom.predict(allom_dir,
                                  dbh = unlist(mw[pnt, c('diam1', 'diam2')]) * cm_per_inch,
                                  pft = unlist(mw[pnt, c('pecan1', 'pecan2')]),
                                  component = allom_component, 
                                  use = 'Bg', # 'mu' in unstable statistically
                                  n = n_allom_samples,
                                  interval = "prediction", 
                                  single.tree = FALSE)
            if(!do_allom_uncertainty) {
                mw[pnt, c('biomass1', 'biomass2')] <- colMeans(pred)  
            } else stop("allometric uncertainty sampling not yet fully coded")                
        } else {
            pred <- allom.predict(allom_dir,
                                  dbh = unlist(mw[pnt, c('diam1')]) * cm_per_inch,
                                  pft = unlist(mw[pnt, c('pecan1')]),
                                  component = allom_component, 
                                  use = 'Bg', # 'mu' in unstable statistically
                                  n = n_allom_samples,
                                  interval = "prediction", 
                                  single.tree = FALSE)
            if(!do_allom_uncertainty) {
                mw[pnt, c('biomass1')] <- mean(pred)
            } else stop("allometric uncertainty sampling not yet fully coded")                
        }
        if(pnt %% 10000 == 0) cat("Finished row ", pnt, "\n")
    }
} else {
    mw <- mw %>% add_cells_to_dataset()
    cells <- unique(mw$cell)
    for(cc in seq_along(cells)) {
        wh <- which(mw$cell == cells[cc])
        two_tree_points <- wh[mw$num_trees[wh] == 2]
        one_tree_points <- wh[mw$num_trees[wh] == 1]
        if(length(two_tree_points)) {
            pred <- allom.predict(allom_dir,
                                  dbh = unlist(mw[two_tree_points, c('diam1', 'diam2')]) * cm_per_inch,
                                  pft = unlist(mw[two_tree_points, c('pecan1', 'pecan2')]),
                                  component = allom_component, 
                                  use = 'Bg', # 'mu' in unstable statistically
                                  n = n_allom_samples,
                                  interval = "prediction", 
                                  single.tree = FALSE)
            if(!do_allom_uncertainty) {
                mw[two_tree_points, c('biomass1', 'biomass2')] <- colMeans(pred)  
            } else stop("allometric uncertainty sampling not yet fully coded")
        }
        if(length(one_tree_points)) {
            pred <- allom.predict(allom_dir,
                                  dbh = unlist(mw[one_tree_points, 'diam1']) * cm_per_inch,
                                  pft = unlist(mw[one_tree_points, 'pecan1']),
                                  component = allom_component, 
                                  use = 'Bg', # 'mu' in unstable statistically
                                  n = n_allom_samples,
                                  interval = "prediction", 
                                  single.tree = FALSE)
            if(!do_allom_uncertainty) {
                mw[one_tree_points, 'biomass1'] <- colMeans(pred)  ## for moment use posterior mean
            } else stop("allometric uncertainty sampling not yet fully coded")                
        }
        if(cc %% 100 == 0) cat("Finished cell index ", cc, "\n")
    }
}
cat("Note: currently using point estimate for individual tree biomass.\n")

## check for outliers in individual biomass draws

## Actually for now do include unknown tree biomass.
## For unknown tree, we could just use some overall average allometry, but perhaps for
## simplicity just don't determine a biomass. There are 398 Unknown tree in tree1 and 370 in tree2.
## This needs to be done because we set a default allometry for 'Unknown tree' above and want
## to undo the effect of that.
if(FALSE) {
    mw <- mw %>% mutate(biomass1 = ifelse(L3_tree1 == "Unknown tree", NA, biomass1),
                        biomass2 = ifelse(L3_tree2 == "Unknown tree", NA, biomass2))
}

## This avoids issue of averaging a 0 biomass when it is really a missing biomass
mw <- mw %>% mutate(biomass1 = ifelse(!is.na(diam1) & diam1 == 0, NA, biomass1),
                    biomass2 = ifelse(!is.na(diam2) & diam2 == 0, NA, biomass2))

## based on issue #39 we are not omitting trees below the veil line
if(FALSE) { 
    mw <- mw %>% mutate(biomass1 = ifelse(!is.na(diam1) & diam1 < diameter_cutoff_inches, NA, biomass1),
                        biomass2 = ifelse(!is.na(diam2) & diam2 < diameter_cutoff_inches, NA, biomass2))
}

if(shared_params_in_cell) {
    save(mw, file = file.path(interim_results_dir, 'point_with_biomass_shared.Rda'))
} else save(mw, file = file.path(interim_results_dir, 'point_with_biomass.Rda'))
