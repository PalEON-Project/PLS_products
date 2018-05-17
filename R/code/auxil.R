library(raster)
base_raster <- raster(xmn = -71000, xmx = 2297000, ncols = 296,
                      ymn = 58000,  ymx = 1498000, nrows = 180,
                      crs = '+init=epsg:3175')
base_raster <- setValues(base_raster, 1:ncell(base_raster))

convert_to_NA <- function(x, missingCodes = c(88888, 99999)) {
    x[x %in% missingCodes] <- NA
    return(x)
}

add_cells_to_dataset <- function(data) {
    library(dplyr)
    points <- data %>% dplyr::select(x,y)
    coordinates(points) <- ~x+y
    proj4string(points) <- CRS('+init=epsg:3175')
    data <- data %>% mutate(cell = raster::extract(base_raster, points))
    return(data)
}

get_angle_inil <- function(bearings, degrees) {
    ##  This function converts the bearing and degree information
    ##  into 360-degree azimuth values

    angle <- matrix(NA, nrow(degrees), ncol(degrees))

    if(any(is.na(bearings)))
        stop("get_angle_inil: Found NA values for bearings and code not set to handle this.")

    card_north <- c('N_NA', 'NA_N', 'N')
    card_east <- c('E_NA', 'NA_E', 'E')
    card_south <- c('S_NA', 'NA_S', 'S')
    card_west <- c('W_NA', 'NA_W', 'W')
    north <- bearings %in% card_north & degrees == 0
    east <- bearings %in% card_east & degrees == 0
    south <- bearings %in% card_south & degrees == 0
    west <- bearings %in% card_west & degrees == 0

    if(any(bearings %in% c(card_north, card_south, card_east, card_west) & degrees != 0))
        stop("get_angle_inil: Found non-zero degrees for cardinal bearings.")
    
    ## set the unidirectional angles to a degrees 360 
    angle[ north ] <- 0
    angle[ south ] <- 180
    angle[ east  ] <- 90 
    angle[ west  ] <- 270

    ne <- bearings == 'N_E'
    se <- bearings == 'S_E'
    sw <- bearings == 'S_W'
    nw <- bearings == 'N_W'

    ## we trust the quadrant info implied by the bearing info and don't shift angles to new quadrants
    degrees[degrees < 0 | degrees > 90] <- 45
    
    ##  convert all angles from a 90-degree value within each quadrat to a 360 degree value.
    ## 'ne' quadrant is angles of 0-90 starting at north, so doesn't need changing.
    angle[ ne ] <- degrees[ ne ]
    ## 'se' quadrant is angles of 0-90 starting at south  
    angle[ se ] <- 180 - degrees[ se ]
    ## 'sw' quadrant is angles of 0-90 starting at south
    angle[ sw ] <- 180 + degrees[ sw ]
    ## 'nw' quadrant is angles of 0-90 starting at north
    angle[ nw ] <- 360 - degrees [ nw ]

    ## should handle cases with 88888 and 99999 values that indicate unknown azimuth
    angle[ !north & !south & !west & !east & !ne & !nw & !se & !sw ] <- NA

    ## a few hundred points with N_W and angle of 0; need these to be 0 so don't create a fifth quadrant with floor(angle/90) later
    angle[ angle == 360 ] <- 0
   return(angle)
}

get_angle_umw <- function(azimuth) {
  #  This function is used in 'step.one.clean.bind_v1.1.R', it converts the
  #  text azimuth strings to numeric, 360 degree values.
  #  This is the vector that will store the values.
   angle <- matrix(as.numeric(NA), nrow(azimuth), ncol(azimuth))
  
  #  This is a special case, generally where the tree is plot center.
    angle[azimuth == '0'] <- 0

    # these will be assigned angle of 45
    azimuth <- gsub("XX", "  ", azimuth)
    azimuth <- gsub("--", "  ", azimuth)
    
  
  #  This is a short function that takes cares of NAs in boolean functions, it's
  #  just a simple wrapper for the boolean function that sets the NA values in
  #  the vector to FALSE.
  fx.na <- function(x) { x[ is.na( x ) ] <- FALSE; x }
  
  #  Given the text azimuths in the dataset, return the quadrant values.
  #  This gives a boolean index of the quadrant

  north <- fx.na(azimuth %in% c('N', 'NORT', 'N0RT'))
  south <- fx.na(azimuth %in% c('S', 'SOUT', 'S0UT'))
  east <- fx.na(azimuth %in% c('E', 'EAST'))
  west <- fx.na(azimuth %in% c('W', 'WEST'))
    
  angle[ north ] <- 0
  angle[ south ] <- 180
  angle[ east  ] <- 90 
  angle[ west  ] <- 270

  card <- c('NORT','EAST','SOUT','WEST','N0RT','S0UT')
  n <- fx.na( regexpr('N', azimuth) > 0 & !azimuth %in% card)
  e  <- fx.na( regexpr('E', azimuth) > 0 & !azimuth %in% card)
  s <- fx.na( regexpr('S', azimuth) > 0 & !azimuth %in% card)
  w <-  fx.na( regexpr('W', azimuth) > 0 & !azimuth %in% card)

    ne <- n & e
    se <- s & e
    sw <- s & w
    nw <- n & w
  #  The cell is in a quadrant, regardless of which.
  quad <- ne | se | sw | nw
  
 
  #  The problem is that some notes have either N04E, N 4E or N4E, or NE!
  strlen <- nchar(azimuth)
  strlen[is.na(azimuth)] <- NA
  
  angle[quad & strlen == 2] <- 45
    angle[quad & strlen == 3] <- as.numeric(substr(azimuth[ quad & strlen == 3 ], 2, 2))

    ## a few O that are presumably 0
    azimuth <- gsub("O", "0", azimuth, fixed = TRUE)

    ## numeric cases; non-numeric such as S6EW will default to NA
    chars <- substr(azimuth, 2, 3)
    wh <- intersect(which(quad & strlen == 4 & chars != '  '), grep("[0-9 ]{2}", chars, perl = TRUE))
    angle[wh] <- as.numeric(chars[wh])
  
  # Special case of double spaces in the azimuth.
    angle[quad & strlen == 4 & chars == '  '] <- 45

    # all of these cases appear safe to simply assume in quadrant
    angle[quad & strlen >= 5] <- 45
    

    ## we trust the quadrant info implied by the bearing info and don't shift angles to new quadrants
    angle[quad & (angle < 0 | angle > 90)] <- 45

    angle[ne] <- angle[ne]
  angle[se] <- 180 - angle[se]
  angle[sw] <- 180 + angle[sw]
  angle[nw] <- 360 - angle[nw]
  
  return(angle)
  
}

calc_stem_density <- function(data, corr_factors, use_phi =  TRUE) {

    corr <- data %>% dplyr::select(state, surveyyear, corner, sectioncorner, point) %>%
        mutate(state = ifelse(grepl("MI", state), "MI", state)) %>%   ## convert soMI, noMI, noMI_extra to MI
        left_join(corr_factors, by = c('state' = 'state', 'surveyyear' = 'year',
                                       'corner' = 'corner',
                                       'sectioncorner' = 'sectioncorner',
                                       'point' = 'point'))
    if(use_phi) {
        corr <- corr %>% dplyr::select(kappa, theta, zeta, phi)
    } else corr <- corr %>% dplyr::select(kappa, theta, zeta)

    data <- data %>% mutate(full_corr = apply(as.matrix(corr), 1, prod))

    ## only determine in same quad if have definitive info to that effect, not if azimuth is missing
    quad1 <- floor(data$az1/90)
    quad2 <- floor(data$az2/90)
    same_quad <- !is.na(quad1) & !is.na(quad2) & quad1 == quad2

    density <- rep(as.numeric(NA), nrow(data))

    ## fix upper radius searched and take 1 tree per that area as density
    max_dist_meters <- max_distance_surveyed * meters_per_chain
    density[data$num_trees == 1] <- (1 / (pi * max_dist_meters^2)) * meters_sq_per_ha
    
    density[data$num_trees == 0] <- 0

    ## We use points where first tree (but not also second tree) has distance of 0
    to_calc <- data$num_trees == 2 & (!(data$dist1 == 0 & data$dist2 == 0))
    sub <- data[to_calc, ]

    diam <- sub$diam1
    diam[is.na(diam)] <- 10
    ## distance in meters: distance plus one-half diameter of tree
    distm1 <- sub$dist1 * meters_per_chain + 0.5 * diam * cm_per_inch / cm_per_m
    diam <- sub$diam2
    diam[is.na(diam)] <- 10
    ## distance in meters: distance plus one-half diameter of tree
    distm2 <- sub$dist2 * meters_per_chain + 0.5 * diam * cm_per_inch / cm_per_m

    ##  From the formula,
    ##  lambda = kappa * theta * (q - 1)/(pi * n) * (q / sum_(1:q)(r^2))
    ##  here, n is equal to 1.
    ##  units are in stems / m^2
    q <- 2

    morisita <- ((q - 1) / pi) * (q / (distm1^2 + distm2^2)) * sub$full_corr 
    morisita <- morisita * meters_sq_per_ha

    density[to_calc] <- morisita

    density[same_quad] <- NA

    density[density > max_density] <- max_density

    return(density)
}
