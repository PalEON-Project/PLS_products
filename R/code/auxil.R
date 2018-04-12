convert_to_NA <- function(x, missingCodes = c(88888, 99999)) {
    x[x %in% missingCodes] <- NA
    return(x)
}

get_angle_inil <- function(bearings, degrees) {
    ##  This function converts the bearing and degree information
    ##  into 360-degree azimuth values

    angle <- degrees

    if(any(is.na(bearings)))
        stop("get_angle_inil: Found NA values for bearings and code not set to handle this.")

    card_north <- c('N_NA', 'NA_N', 'N')
    card_south <- c('E_NA', 'NA_E', 'E')
    card_east <- c('S_NA', 'NA_S', 'S')
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

    ##  convert all angles from a 90-degree value within each quadrat to a 360 degree value.
    ## 'ne' quadrant is angles of 0-90 starting at north, so doesn't need changing.
    ## 'se' quadrant is angles of 0-90 starting at south  
    angle[ se ] <- 180 - degrees[ se ]
    ## 'sw' quadrant is angles of 0-90 starting at south
    angle[ sw ] <- 180 + degrees[ sw ]
    ## 'nw' quadrant is angles of 0-90 starting at north
    angle[ nw ] <- 360 - degrees [ nw ]

    ## should handle cases with 88888 and 99999 values that indicate unknown azimuth
    angle[ !north & !south & !west & !east & !ne & !nw & !se & !sw ] <- NA

    ## only a few of these but they would put the azimuth in a different quadrant
    ## TODO: if we decide to trust the bearing but not degrees, just assign 45,135,225,315 
    angle[degrees >= 90] <- NA

   return(angle)
  
}

get_angle_umw <- function(azimuth) {
  #  This function is used in 'step.one.clean.bind_v1.1.R', it converts the
  #  text azimuth strings to numeric, 360 degree values.
  #  This is the vector that will store the values.
   angle <- matrix(as.numeric(NA), nrow(azimuth), ncol(azimuth))
  
  #  This is a special case, generally where the tree is plot center.
  angle[azimuth == '0'] <- 0
  
  #  This is a short function that takes cares of NAs in boolean functions, it's
  #  just a simple wrapper for the boolean function that sets the NA values in
  #  the vector to FALSE.
  fx.na <- function(x) { x[ is.na( x ) ] <- FALSE; x }
  
  #  Given the text azimuths in the dataset, return the quadrant values.
  #  This gives a boolean index of the quadrant

  north <- fx.na(azimuth %in% c('N', 'NORT'))
  south <- fx.na(azimuth %in% c('S', 'SOUTH'))
  east <- fx.na(azimuth %in% c('E', 'EAST'))
  west <- fx.na(azimuth %in% c('W', 'WEST'))
    
  angle[ north ] <- 0
  angle[ south ] <- 180
  angle[ east  ] <- 90 
  angle[ west  ] <- 270
    
  n <- fx.na( regexpr('N', azimuth) > 0 & azmith != 'NORT')
  e  <- fx.na( regexpr('E', azimuth) > 0 & azimuth != 'EAST')
  s <- fx.na( regexpr('S', azimuth) > 0 & azimuth != 'SOUT')
  w <-  fx.na( regexpr('W', azimuth) > 0 & azimuth != 'WEST')

  #  The cell is in a quadrant, regardless of which.
  quad <- ( (n&e) | (s&e) | (s&w) | (n&w) )
  
 
  #  The problem is that some notes have either N04E, N 4E or N4E, or NE!
  strlen <- nchar(azimuth)
  strlen[is.na(azimuth)] <- NA
  
  angle[quad & strlen == 2] <- 45
  angle[quad & strlen == 3] <- as.numeric(substr(azimuth[ quad & strlen == 3 ], 2, 2))
  angle[quad & strlen == 4 & !substr(azimuth, 2, 3) == '  '] <- 
    as.numeric(substr(azimuth[quad & strlen == 4 & !substr(azimuth, 2, 3) == '  '], 2, 3))
  
  # Special case of double spaces in the azimuth.
    angle[quad & strlen == 4 & substr(azimuth, 2, 3) == '  '] <- 45

    ## TEMP, while we get feedback from Charlie/Simon
    angle[angle >= 90] <- NA

    ## cases such as S--W, S69.3W, 1991, N9999W will default to NA
    ## TODO: check back to see if we want to treat these as known quadrants
  angle[ se ] <- 180 - angle[ se ]
  angle[ sw ] <- 180 + angle[ sw ]
  angle[ nw ] <- 360 - angle [ nw ]
  
  return(angle)
  
}
