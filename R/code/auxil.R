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

