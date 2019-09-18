## Estimate tree-level basal area using diameter of trees 

library(dplyr)
library(assertthat)

load(file.path(interim_results_dir, 'point_with_density.Rda'))

## calculate basal area in m^2
calc_basalarea <- function(diam_inch)
    return( pi * ((diam_inch * cm_per_inch / cm_per_m) / 2)^2 )

mw <- mw %>% mutate(basalarea1 = calc_basalarea(diam1), basalarea2 = calc_basalarea(diam2))
    
max_basalarea <- pi * ((150 * cm_per_inch / cm_per_m) / 2)^2

assert_that(min(mw$basalarea1, na.rm = TRUE) >= 0 & min(mw$basalarea2, na.rm = TRUE) >= 0 &
            max(mw$basalarea1, na.rm = TRUE) < max_basalarea & max(mw$basalarea2, na.rm = TRUE) < max_basalarea,
            msg = "extreme basal area values")

## This avoids issue of averaging a 0 basal area when it is really a missing basal area
mw <- mw %>% mutate(basalarea1 = ifelse(!is.na(diam1) & diam1 == 0, NA, basalarea1),
                    basalarea2 = ifelse(!is.na(diam2) & diam2 == 0, NA, basalarea2))

save(mw, file = file.path(interim_results_dir, 'point_with_basalarea.Rda'))
