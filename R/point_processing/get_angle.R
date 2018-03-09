#  This function deconvolves the azimuth as stored in the 

get_angle <- function(azimuth) {
  #  This function is used in 'step.one.clean.bind_v1.1.R', it converts the
  #  text azimuth strings to numeric, 360 degree values.
  #  This is the vector that will store the values.
  angl <- rep(NA, length(azimuth))
  
  #  This is a special case, generally where the tree is plot center.
  angl[azimuth == '0'] <- 0
  
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
  
  angl[ uni & north ] <- 0
  angl[ uni & south ] <- 180
  angl[ uni & east  ] <- 90 
  angl[ uni & west  ] <- 270
  
  #  The problem is that some notes have either N04E, N 4E or N4E, or NE!
  strlen <- nchar(azimuth)
  strlen[is.na(azimuth)] <- NA
  
  angl[quad & strlen == 2] <- 45
  angl[quad & strlen == 3] <- as.numeric(substr(azimuth[ quad & strlen == 3 ], 2, 2))
  angl[quad & strlen == 4 & !substr(azimuth, 2, 3) == '  '] <- 
    as.numeric(substr(azimuth[quad & strlen == 4 & !substr(azimuth, 2, 3) == '  '], 2, 3))
  
  # Special case of double spaces in the azimuth.
  angl[quad & strlen == 4 & substr(azimuth, 2, 3) == '  '] <- 45
  
  #  Another set of special cases:
  angl[ fx.na(azimuth == 'NORT') ] <- 0
  angl[ fx.na(azimuth == 'EAST') ] <- 90
  angl[ fx.na(azimuth == 'WEST') ] <- 270
  angl[ fx.na(azimuth == 'SOUT') ] <- 180
  
  angl[ se ] <- 180 - angl[ se ]
  angl[ sw ] <- 180 + angl[ sw ]
  angl[ nw ] <- 360 - angl [ nw ]
  
  return(angl)
  
}
