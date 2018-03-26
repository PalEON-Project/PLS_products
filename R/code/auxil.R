convert_to_NA <- function(x, missingCodes = c(88888, 99999)) {
    x[x %in% missingCodes] <- NA
    return(x)
}

get_angle_inil <- function(bearings, degrees, dists) {
  ##  This function converts the
  ##  text azimuth strings to numeric, 360 degree values.
  ##  This is the vector that will store the values.
  angle <- degrees
  
  ##  This is a special case, generally where the tree is plot center.
    ## CJP: this doesn't make any sense - angle is equal to degrees, so
    ## if degrees is 0 then angle already is 0
  if(FALSE) {
      angle[degrees == 0 & dists == 0] <- 0
  }
    
  ##  This sets NAs to FALSE
  fx.na <- function(x) { x[ is.na( x ) ] <- FALSE; x }
  
  ##  Given the text azimuths in the dataset, return the quadrant values.
  ##  This gives a boolean index of the quadrant
  
  north <- fx.na(bearings == 'N_NA'| bearings == 'NA_N' | bearings =='N') #|bearings =='N99999'|bearings =='N88888')
  east <- fx.na(bearings == 'NA_E' | bearings =="E_NA" | bearings =='E') #|bearings =='99999E'|bearings =='88888E')
  south <- fx.na(bearings == 'S_NA' | bearings =="NA_S"| bearings == 'S') #|bearings =='S99999'|bearings =='S88888')
  west <- fx.na(bearings == 'NA_W' | bearings =="W_NA" | bearings == 'W') #|bearings =='99999W'|bearings =='88888W')
  
  ##north <- fx.na( regexpr('N', bearings) > 0 )
  ##east  <- fx.na( regexpr('E', bearings) > 0 | bearings == 'EAST')
  ##south <- fx.na( regexpr('S', bearings) > 0 | bearings == 'SOUTH')
  ##west <-  fx.na( regexpr('W', bearings) > 0 | bearings == 'WEST')

    ## CJP: this doesn't make any sense - by how north,east,south,west
    ## defined above, the first condition in these four lines will never be TRUE
  if(FALSE) {
      ne <- fx.na( (north & east) | bearings == 'N_E')
      se <- fx.na( (south & east) | bearings == 'S_E')
      sw <- fx.na( (south & west) | bearings == 'S_W')
      nw <- fx.na( (north & west) | bearings == 'N_W') 
  }
    
      ne <- fx.na(bearings == 'N_E')
      se <- fx.na(bearings == 'S_E')
      sw <- fx.na(bearings == 'S_W')
      nw <- fx.na(bearings == 'N_W') 

    ##  The cell is in a quadrant, regardless of which.
  quad <- ne | se | sw | nw
  
  ##  Special case of the trees with a unidirectional direction.
    uni  <- (!quad) & !(north | south | east | west)
    ## CJP: what does unidirectional mean - based on the code here, aren't these cases where we have
    ## NA for both bearing and bearingdir in the original data file?

    ## CJP: this doesn't make any sense - uni is only TRUE if none of north/south/east/west are,
    ## so none of these condition are ever true
    if(FALSE) {
        angle[ uni & north ] <- 0
        angle[ uni & south ] <- 180
        angle[ uni & east  ] <- 90 
        angle[ uni & west  ] <- 270
    }
    ## CJP: as a result, quad,ne,se,nw,sw do not need to be computed
  
   ##  Another set of special cases:
 
    angle[ north ] <- degrees[ north ] ## CJP: this has no effect since angle started out equal to degrees
    ## CJP: would be good to have comment that explains what this is doing.
  angle[ east ] <- 180 - degrees[ east ]
  angle[ south ] <- 180 + degrees[ south ]
  angle[ west ] <- 360 - degrees [ west ]
  
  return(angle)
  
}

