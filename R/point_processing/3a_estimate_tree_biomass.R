## using code from CJP fia_products estimate biomass using PEcAn at tree level

## get allometry samples via PEcAn

library(dplyr)
library(readr)
library(PEcAn.allometry)
library(PEcAn.logger)

cols_conv <- cols(
  level3a = col_character(),
  pecan_allometry_spcd = col_character(),  ## multiple spcd separated by ; for non-exact matches
  pecan_allometry_common_names = col_character()
)


load('point_with_density.Rda')

## why is MASS::select being loaded (why MASS -- comes from Pecan)

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

## updated more lumped conversions based on initial results commented about below
## TODO: this should be done in level3a_to_fia
to_convert <- c('762','970')
mw <- mw %>% mutate(pecan1 = ifelse(pecan1 %in% to_convert, '318,802,541,731', pecan1),
                    pecan2 = ifelse(pecan2 %in% to_convert, '318,802,541,731', pecan2))
to_convert <- c('491', '701')
mw <- mw %>% mutate(pecan1 = ifelse(pecan1 %in% to_convert, '315,471,491,319,355,356,357,931,935,391,701,660,761,500,501,502,763,662,663,765,766', pecan1),
                    pecan2 = ifelse(pecan2 %in% to_convert, '315,471,491,319,355,356,357,931,935,391,701,660,761,500,501,502,763,662,663,765,766', pecan2))
to_convert <- c('540')
mw <- mw %>% mutate(pecan1 = ifelse(pecan1 %in% to_convert, '541,543,544', pecan1),
                    pecan2 = ifelse(pecan2 %in% to_convert, '541,543,544', pecan2))


## placeholder so that NA values don't cause problems in allometry fitting
mw <- mw %>% mutate(pecan1 = ifelse(is.na(pecan1), '318,802,541,731', pecan1),
                    pecan2 = ifelse(is.na(pecan2), '318,802,541,731', pecan2))

## Fit allometry models for all taxa

# actually we should be fitting for the pecan_allometry_spcd in the data not in the conversion table but probably fine

unique_pecan_allom <- unique(c(mw$pecan1, mw$pecan2))
unique_pecan_allom <- unique_pecan_allom[!is.na(unique_pecan_allom)]  ## from Unknown tree
pecan_taxa <- lapply(seq_len(length(unique_pecan_allom)), function(i) 
    data.frame(spcd = as.numeric(strsplit(unique_pecan_allom[i], split = ',')[[1]])))
names(pecan_taxa) <- unique_pecan_allom

## first check which allometries are already in the allom_dir and don't redo
## TBD
allom_stats = load.allom(allom_dir)

##allom_stats = AllomAve(pecan_taxa, ngibbs=1000, components = 6, outdir = allom_dir, dmin = 10, dmax = 150)

## create via loop rather than passing list of dataframes because of PEcAn bug -
## errors out for spcd=701
for(i in seq_along(pecan_taxa)) {
    if(!names(pecan_taxa)[i] %in% names(allom_stats)) {
        cat("Fitting ", names(pecan_taxa)[i], ".\n")
        allom_stats[[names(pecan_taxa)[i]]] <- try(AllomAve(pecan_taxa[i], ngibbs=1000, components = 6,
                             outdir = allom_dir, dmin = 10, dmax = 150))
}}

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


## site effects seem to product crazy allometries - mu0 and mu1 have outliers and taus can be big; how many allometries do we need for stability?

## We want to run the allometry prediction for all trees of a given taxon at a given point simultaneously because this allows for shared allometric parameters (but different tree effects) for trees of the same taxon. In this case we simply run the prediction for all trees at a point at once.
## arguably we should use shared parameters for all trees in a given cell, but hard to know what spatial scale to assume homogeneity at

mw$biomass1 <- mw$biomass2 <- rep(NA, nrow(mw))

## This takes ?? (3.5) hours on smeagol on one core.
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
            mw[pnt, c('biomass1', 'biomass2')] <- colMeans(pred)  ## for moment use posterior mean
        } else {
            pred <- allom.predict(allom_dir,
                                  dbh = unlist(mw[pnt, c('diam1')]) * cm_per_inch,
                                  pft = unlist(mw[pnt, c('pecan1')]),
                                  component = allom_component, 
                                  use = 'Bg', # 'mu' in unstable statistically
                                  n = n_allom_samples,
                                  interval = "prediction", 
                                  single.tree = FALSE)
            mw[pnt, c('biomass1')] <- mean(pred)
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
            mw[two_tree_points, c('biomass1', 'biomass2')] <- colMeans(pred)  ## for moment use posterior mean
        }
        if(length(one_tree_points)) {
            pred <- allom.predict(allom_dir,
                                  dbh = unlist(mw[one_tree_points, c('diam1', 'diam2')]) * cm_per_inch,
                                  pft = unlist(mw[one_tree_points, c('pecan1', 'pecan2')]),
                                  component = allom_component, 
                                  use = 'Bg', # 'mu' in unstable statistically
                                  n = n_allom_samples,
                                  interval = "prediction", 
                                  single.tree = FALSE)
            mw[one_tree_points, c('biomass1', 'biomass2')] <- colMeans(pred)  ## for moment use posterior mean
        }
        if(cc %% 100 == 0) cat("Finished cell index ", cc, "\n")
    }
}

## check for outliers in individual biomass draws


## For unknown tree, we could just use some overall average allometry, but perhaps for
## simplicity just don't determine a biomass. There are 398 Unknown tree in tree1 and 370 in tree2.
## This needs to be done because we set a default allometry for 'Unknown tree' above and want
## to 'undo' the effect of that.
mw <- mw %>% mutate(biomass1 = ifelse(L3_tree1 == "Unknown tree", NA, biomass1),
                    biomass2 = ifelse(L3_tree2 == "Unknown tree", NA, biomass2))


mw <- mw %>% mutate(biomass1 = ifelse(!is.na(diam1) & diam1 < diameter_cutoff_inches, NA, biomass1),
                    biomass2 = ifelse(!is.na(diam2) & diam2 < diameter_cutoff_inches, NA, biomass2))

mn <- mw %>% dplyr::select(biomass1, biomass2) %>% rowMeans(na.rm = TRUE)

mw <- mw %>% mutate(biomass = ifelse(num_trees == 0, 0,
                              ifelse(num_trees == 1, biomass1, mn)))
mw <- mw %>% mutate(biomass = ifelse(is.nan(biomass), NA, biomass))


save(mw, file = 'point_with_biomass.Rda')
