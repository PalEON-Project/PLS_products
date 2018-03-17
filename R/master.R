## Full workflow for PLS products for PalEON midwest region: biomass, density, basal area
## Chris Paciorek
## 2018-02-15

## note required input files (download from PalEON Wiki)


## upper_midwest_minnesota_file="umw_mn_pointdata_v1.0.csv"
## upper_midwest_wisconsin_file="umw_wi_pointdata_v1.0.csv"
## upper_midwest_michigan_file="umw_mi_pointdata_v1.0.csv"
## southern_michigan_file="southernmi_projected_v1.2.csv"
## indiana_file="ndinpls_v1.7.csv"
## illinois_file="ndilpls_v1.8-2.csv"
## taxa_conversion_file="level0_to_level3a_v0.4-7.csv"

## setup R packages (with version control) 

checkpoint::checkpoint("2018-01-18")

library(devtools)
install_github("PecanProject/pecan",subdir="base/logger", ref = pecan_base_logger_commit)
install_github("PecanProject/pecan",subdir="modules/allometry", ref = pecan_modules_allometry_commit)

## get configuration variables

source('config')

## setup directories and basic data

states = c('MN', 'WI', 'MI', 'IL', 'IN', 'OH', 'PA', 'NY', 'NJ', 'ME', 'VT', 'MA', 'CT', 'NH', 'RI')

paleon_regions_west <- c(2,3,5,6,11,12) ## MI/IN and west
paleon_regions <- 1:12

paleon_states_west <- c(17, 18, 26, 27, 55) 

# check this:
# excluded_level3s_OH <- c('Dogwood') # only 20 trees >= 20 cm
## Note that all Atlantic white cedar is mistakenly classified as cedar/juniper; have asked Jody to change this

if(raw_data_dir == "")
    raw_data_dir <- file.path("..", "data", "points") 
if(conversions_data_dir == "")
    conversions_data_dir <- file.path("..", "data", "conversions") 
if(allom_dir == "")
    allom_dir <- file.path("..", "data", "allom")
if(code_dir == "")
    code_dir <- file.path("code")

## source files with R functions

code_files <- list.files(code_dir, full.names = TRUE)
sapply(code_files, source)

# TMP: for seeing all cols instead of tibble
adf <- as.data.frame
