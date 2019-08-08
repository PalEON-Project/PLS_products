## Estimate tree-level biomass using diameter of trees and PEcAn allometries

## Run-time: approximately 60 minutes when using allometry parameters shared by cell with 25 samples.

library(dplyr)
library(readr)
library(PEcAn.allometry)
library(PEcAn.logger)
library(assertthat)

load(file.path(interim_results_dir, 'point_with_density.Rda'))

## Read in conversion info to go from FIA taxa to PalEON-determined aggregations of taxa
## for which we can estimate statistically stable allometric relationships.

cols_conv <- cols(
  level3a = col_character(),
  pecan_allometry_spcd = col_character(),  ## multiple spcd separated by ; for non-exact matches
  pecan_allometry_common_names = col_character()
)

## MASS masks select, so need dplyr::select (MASS -- comes from Pecan)

if(use_agb) {

    if(!file.exists(file.path(conversions_data_dir, pls_to_chojnacky_conversion_file))) 
        unzip(file.path(conversions_data_dir, pls_to_chojnacky_zipfile), exdir = conversions_data_dir)

    taxa_conversion <- read_csv(file.path(conversions_data_dir, pls_to_chojnacky_conversion_file),
                                col_types = cols_conv) %>%
        dplyr::select(level3a, beta0, beta1)

    mw <- mw %>% left_join(taxa_conversion, by = c("L3_tree1" = "level3a")) %>%
        rename(int1 = beta0, slope1 = beta1) %>%
        left_join(taxa_conversion, by = c("L3_tree2" = "level3a")) %>%
        rename(int2 = beta0, slope2 = beta1)

    predict_biomass <- function(dbh_inch, b0, b1) {
        return(exp(b0+b1*log(dbh_inch * cm_per_inch)))
    }

    mw <- mw %>% mutate(biomass1 = predict_biomass(diam1, int1, slope1)) %>% 
        mutate(biomass2 = predict_biomass(diam2, int2, slope2))


}  else {
    
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


## Treat Unknown tree as generic hardwood.
## Even if we want to omit Unknown tree, we need a placeholder so that NA
## (from Unknown tree) values don't cause problems in allometry fitting
assert_that(sum(is.na(mw$pecan1) & !is.na(mw$L3_tree1) & mw$L3_tree1 != "No tree") ==
            sum(mw$L3_tree1 == "Unknown tree", na.rm = TRUE),
            msg = "Not all unknown pecan1 values are unknown tree")
assert_that(sum(is.na(mw$pecan2) & !is.na(mw$L3_tree2) & mw$L3_tree2 != "No tree") ==
            sum(mw$L3_tree2 == "Unknown tree", na.rm = TRUE),
            msg = "Not all unknown pecan1 values are unknown tree")
mw <- mw %>% mutate(pecan1 = ifelse(is.na(pecan1), '318,802,541,731', pecan1),
                    pecan2 = ifelse(is.na(pecan2), '318,802,541,731', pecan2))

## Fit allometry models for all taxa
unique_pecan_allom <- unique(c(mw$pecan1, mw$pecan2))
pecan_taxa <- lapply(seq_len(length(unique_pecan_allom)), function(i) 
    data.frame(spcd = as.numeric(strsplit(unique_pecan_allom[i], split = ',')[[1]])))
names(pecan_taxa) <- unique_pecan_allom

## First check which allometries are already in the allom_dir and don't redo
allom_stats = load.allom(allom_dir)

## Create via loop rather than passing list of dataframes because of PEcAn bug 
for(i in seq_along(pecan_taxa)) {
    if(!names(pecan_taxa)[i] %in% names(allom_stats)) {
        cat("Fitting ", names(pecan_taxa)[i], ".\n")
        allom_stats[[names(pecan_taxa)[i]]] <- try(AllomAve(pecan_taxa[i], ngibbs = 1000,
                                                            components = allom_component,
                                                            outdir = allom_dir, dmin = dmin_allom_fit, dmax = dmax_allom_fit))
    }}



## We want to run the allometry prediction for all trees of a given taxon at a given point simultaneously
## because this allows for shared allometric parameters (but different tree effects) for trees of the same taxon.
## In this case we simply run the prediction for all trees at a point at once.

## Actually, arguably we should use shared parameters for all trees in a given cell,
## but hard to know what spatial scale to assume homogeneity at.
## Given the increased computational time of running separately for each point,
## we are currently (based on shared_params_in_cell) doing for all trees in a cell.

mw <- mw %>% mutate(biomass1 = NA, biomass2 = NA)

## This takes ?? (3.5) hours on one core.
## Parallelizing would require data output manipulation because we write into rows of 'mw'
## Note that every call loads in the allometry material for all of the fitted allometries,
## so it would be easy to make this more efficient.
set.seed(1)

## NAs in diam become in NAs in biomass (with warnings)

## Inclusion of site effects seem to product crazy allometries:
## mu0 and mu1 have outliers and taus can be big, so have 'use' be 'Bg'
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
            pred[pred > biomass_max_kg] <- biomass_max_kg
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
            pred[pred > biomass_max_kg] <- biomass_max_kg
            if(!do_allom_uncertainty) {
                mw[pnt, c('biomass1')] <- mean(pred)
            } else stop("allometric uncertainty sampling not yet fully coded")                
        }
        if(pnt %% 10000 == 0) cat("Finished row ", pnt, "\n")
    }
} else {
    mw <- mw %>% add_paleon_grid()
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
            pred[pred > biomass_max_kg] <- biomass_max_kg
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
            pred[pred > biomass_max_kg] <- biomass_max_kg
            if(!do_allom_uncertainty) {
                mw[one_tree_points, 'biomass1'] <- colMeans(pred)  ## for moment use posterior mean
            } else stop("allometric uncertainty sampling not yet fully coded")                
        }
        if(cc %% 100 == 0) cat("Finished cell index ", cc, "\n")
    }
}
cat("Note: currently using point estimate for individual tree biomass.\n")

}
    
assert_that(min(mw$biomass1, na.rm = TRUE) >= 0 & min(mw$biomass2, na.rm = TRUE) >= 0 &
            max(mw$biomass1, na.rm = TRUE) < 1e6 & max(mw$biomass2, na.rm = TRUE) < 1e6,
            msg = "extreme biomass values")

## We considered setting biomass for 'Unknown tree' to NA, but for now we do estimate biomass
## by using generic hardwood allometry (approximately 80% of trees in database are hardwoods).
## There are 388 Unknown tree in tree1 and 362 in tree2.
if(FALSE) {
    ## This needs to be done because we set a default allometry for 'Unknown tree' above and want
    ## to undo the effect of that.
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

fn <- 'point_with_biomass'
if(shared_params_in_cell)
    fn <- paste0(fn, '_shared')
if(use_agb)
    fn <- paste0(fn, '_agb')
fn <-  paste0(fn, '.Rda')

save(mw, file = file.path(interim_results_dir, fn))
