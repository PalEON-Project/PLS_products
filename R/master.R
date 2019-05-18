## Full workflow for PLS products for PalEON midwest region: biomass, density, basal area (eventually)
## Chris Paciorek
## 2018-02-15

## get configuration variables
source('config')

## set up R packages (with version control)

if(use_original_R_package_versions) {  
    ## bug in checkpoint: errors out with 'arguments imply differing number of rows':
    ## work-around is to avoid wget as download method
    if(options('download.file.method') == 'wget')
        options('download.file.method' = NULL)
    checkpoint::checkpoint(checkpoint_date, R.version = R_version)
    ## devtools won't reinstall if not needed
    devtools::install_github("PecanProject/pecan",subdir="base/logger", ref = pecan_base_logger_commit)
    devtools::install_github("PecanProject/pecan",subdir="modules/allometry", ref = pecan_modules_allometry_commit)
}

## These packages are needed by PEcAn.allometry so fake-loading here to force checkpoint to install them (checkpoint scans but doesn't run code...).

if(FALSE) {
    library(MCMCpack)
    library(mvtnorm)
    library(coda)
    library(XML)
}


## geographic subsetting

states = c('MN', 'WI', 'MI', 'IL', 'IN', 'OH', 'PA', 'NY', 'NJ', 'ME', 'VT', 'MA', 'CT', 'NH', 'RI')

## PalEON-defined regions: some PalEON input files are stratified
## by our own regions that cross state boundaries
paleon_regions_west <- c(2,3,5,6,11,12) ## MI/IN and west
paleon_regions <- 1:12

## State FIPS codes:
paleon_states_west <- c(17, 18, 26, 27, 55) 

## setup directories relative to current working directory,
## if directories not set in config file
if(raw_data_dir == "")
    raw_data_dir <- file.path("..", "data", "points") 
if(conversions_data_dir == "")
    conversions_data_dir <- file.path("..", "data", "conversions") 
if(allom_dir == "")
    allom_dir <- file.path("..", "data", "allom")
if(code_dir == "")
    code_dir <- file.path("code")
if(output_dir == "")
    output_dir <- file.path("..", "output")
if(interim_results_dir == "")
    interim_results_dir <- file.path("..", "data", "interim")

## source files with R functions
code_files <- list.files(code_dir, pattern = ".R$", full.names = TRUE)
sapply(code_files, source)

excluded_level3s <- c('Chestnut')
## very few trees so can't model separately

## values used in cross-validation
## caution: k values of 3000,3500 for 'occ' can take a very long time to fit

k_occ_cv <- c(100,250,500,1000,1500,2000,2500, 3000,3500)
## warning("Not using k_occ of 3000 or 3500 at the moment for CV")
k_pot_cv = c(100,250,500,1000,1500,2000,2500,3000,3500)

## The following documents the workflow to do the analyses
## but would not generally be run all at once,
## so embedded in conditional statements.

## One needs to source master.R before running any of these.

## for data download and processing:
if(FALSE) {
    source(file.path("point_processing", "1_clean_merge_point_data.R"))
    source(file.path("point_processing", "2_estimate_point_density.R"))
    source(file.path("point_processing", "3_estimate_tree_biomass.R"))
}

## for biomass modeling
if(FALSE) {
    source(file.path("stat_modeling", "1_setup_biomass.R"))
    ## to determine optimal k_occ, k_pot, run these two lines:
    source(file.path("stat_modeling", "2_cv_total_biomass.R"))
    source(file.path("stat_modeling", "2_cv_taxon_biomass.R"))
    ## based on results, set k_occ_taxon, k_pot_taxon, k_occ_total, k_pot_total in 'config'
    ## then run:
    source(file.path("stat_modeling", "3_fit_total_biomass.R"))
    source(file.path("stat_modeling", "3_fit_taxon_biomass.R"))
    source(file.path("stat_modeling", "4_output_biomass.R"))
}

## for density modeling
if(FALSE) {
    source(file.path("stat_modeling", "1_setup_density.R"))
    ## to determine optimal k_occ, k_pot, run these two lines:
    source(file.path("stat_modeling", "2_cv_total_density.R"))
    source(file.path("stat_modeling", "2_cv_taxon_density.R"))
    ## based on results, set k_occ_taxon, k_pot_taxon, k_occ_total, k_pot_total in 'config'
    ## then run:
    source(file.path("stat_modeling", "3_fit_total_density.R"))
    source(file.path("stat_modeling", "3_fit_taxon_density.R"))
    source(file.path("stat_modeling", "4_output_density.R"))
}

# TMP: for seeing all cols instead of tibble
adf <- as.data.frame
