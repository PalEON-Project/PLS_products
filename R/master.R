## Full workflow for PLS products for PalEON midwest region: biomass, density (not complete here), basal area (eventually)
## Chris Paciorek
## 2018-02-15

## setup R packages (with version control) 

if(FALSE) {  ## turn off during initial development  
    checkpoint::checkpoint("2018-01-18")
    library(devtools)
    install_github("PecanProject/pecan",subdir="base/logger", ref = pecan_base_logger_commit)
    install_github("PecanProject/pecan",subdir="modules/allometry", ref = pecan_modules_allometry_commit)
}

## note required input files in 'config' (download from PalEON Wiki)
## get configuration variables
source('config')

## region/state info

states = c('MN', 'WI', 'MI', 'IL', 'IN', 'OH', 'PA', 'NY', 'NJ', 'ME', 'VT', 'MA', 'CT', 'NH', 'RI')

paleon_regions_west <- c(2,3,5,6,11,12) ## MI/IN and west
paleon_regions <- 1:12

paleon_states_west <- c(17, 18, 26, 27, 55) 

## setup directories relative to current working directory

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

## values used in cross-validation
k_occ_cv <- c(100,250,500,1000,1500,2000,2500)
k_pot_cv = c(100,250,500,1000,1500,2000,2500,3000,3500)


## source files with R functions
code_files <- list.files(code_dir, pattern = ".R$", full.names = TRUE)
sapply(code_files, source)

# TMP: for seeing all cols instead of tibble
adf <- as.data.frame
