convert_to_NA <- function(x, missingCodes = c(88888, 99999)) {
    x[x %in% missingCodes] <- NA
    return(x)
}

get_angle_inil <- function(bearings, degrees, dists) {
    ##  This function converts the
    ##  text azimuth strings to numeric, 360 degree values.
    ##  This is the vector that will store the values.
    angle <- degrees
    
    ##  This sets NAs to FALSE
    fx.na <- function(x) { x[ is.na( x ) ] <- FALSE; x }
    
    ##  Given the text azimuths in the dataset, return the quadrant values.
    ##  This gives a boolean index of the quadrant
    
    north <- fx.na(bearings == 'N_NA'| bearings == 'NA_N' | bearings =='N') 
    east <- fx.na(bearings == 'NA_E' | bearings =="E_NA" | bearings =='E') 
    south <- fx.na(bearings == 'S_NA' | bearings =="NA_S"| bearings == 'S') 
    west <- fx.na(bearings == 'NA_W' | bearings =="W_NA" | bearings == 'W') 
    
    ne <- fx.na(bearings == 'N_E')
    se <- fx.na(bearings == 'S_E')
    sw <- fx.na(bearings == 'S_W')
    nw <- fx.na(bearings == 'N_W') 
    
    ## set the unidirectional angles to a degrees 360 
    angle[ north ] <- 360 # KAH changed this from 0 to 360 to match convention with southwest michigan, where an angle == 0 is a missing/erroreous value
    angle[ south ] <- 180
    angle[ east  ] <- 90 
    angle[ west  ] <- 270

    ##  convert all angles from a 90-degree value within each quadrat to a 360 degree value.
    ## 'ne' quadrant is angles of 0-90 starting at north, so doesn't need changing.
    ## 'se' quadrant is angles of 0-90 starting at south  
    angle[ se ] <- 180 - degrees[ se ]
    ## 'sw' quadrant is angles of 0-90 starting at south
    angle[ sw ] <- 180 + degrees[ sw ]
    ## 'nw' quadrant is angles of 0-90 starting at north
    angle[ nw ] <- 360 - degrees [ nw ]

   return(angle)
  
}

get_angle_umw <- function(azimuth) {
  #  This function is used in 'step.one.clean.bind_v1.1.R', it converts the
  #  text azimuth strings to numeric, 360 degree values.
  #  This is the vector that will store the values.
  angle <- rep(NA, length(azimuth))
  
  #  This is a special case, generally where the tree is plot center.
  angle[azimuth == '0'] <- 0
  
  #  This is a short function that takes cares of NAs in boolean functions, it's
  #  just a simple wrapper for the boolean function that sets the NA values in
  #  the vector to FALSE.
  fx.na <- function(x) { x[ is.na( x ) ] <- FALSE; x }
  
  #  Given the text azimuths in the dataset, return the quadrant values.
  #  This gives a boolean index of the quadrant
  
  north <- fx.na( regexpr('N', azimuth) > 0 )
  east  <- fx.na( regexpr('E', azimuth) > 0 | azimuth == 'EAST')
  south <- fx.na( regexpr('S', azimuth) > 0 | azimuth == 'SOUT')
  west <-  fx.na( regexpr('W', azimuth) > 0 | azimuth == 'WEST')
  
  ne <- fx.na( (north & east) | azimuth == 'N  E')
  se <- fx.na( (south & east) | azimuth == 'S  E')
  sw <- fx.na( (south & west) | azimuth == 'S  W')
  nw <- fx.na( (north & west) | azimuth == 'N  W') 
  
  #  The cell is in a quadrant, regardless of which.
  quad <- ne | se | sw | nw
  
  #  Special case of the trees with a unidirectional direction.
  uni  <- (!quad) & (north | south | east | west) & (nchar(azimuth) == 1)
  
  angle[ uni & north ] <- 360
  angle[ uni & south ] <- 180
  angle[ uni & east  ] <- 90 
  angle[ uni & west  ] <- 270
  
  #  The problem is that some notes have either N04E, N 4E or N4E, or NE!
  strlen <- nchar(azimuth)
  strlen[is.na(azimuth)] <- NA
  
  angle[quad & strlen == 2] <- 45
  angle[quad & strlen == 3] <- as.numeric(substr(azimuth[ quad & strlen == 3 ], 2, 2))
  angle[quad & strlen == 4 & !substr(azimuth, 2, 3) == '  '] <- 
    as.numeric(substr(azimuth[quad & strlen == 4 & !substr(azimuth, 2, 3) == '  '], 2, 3))
  
  # Special case of double spaces in the azimuth.
  angle[quad & strlen == 4 & substr(azimuth, 2, 3) == '  '] <- 45

    ## CJP: why is this not handled by code such as "azimuth=='EAST'" above? check on this
  #  Another set of special cases:
  angle[ fx.na(azimuth == 'NORT') ] <- 360
  angle[ fx.na(azimuth == 'EAST') ] <- 90
  angle[ fx.na(azimuth == 'WEST') ] <- 270
  angle[ fx.na(azimuth == 'SOUT') ] <- 180
  
  angle[ se ] <- 180 - angle[ se ]
  angle[ sw ] <- 180 + angle[ sw ]
  angle[ nw ] <- 360 - angle [ nw ]
  
  return(angle)
  
}
