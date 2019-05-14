library(raster)
base_raster <- raster(xmn = -71000, xmx = 2297000, ncols = 296,
                      ymn = 58000,  ymx = 1498000, nrows = 180,
                      crs = '+init=epsg:3175')
base_raster <- setValues(base_raster, 1:ncell(base_raster))

convert_to_NA <- function(x, missingCodes = c(88888, 99999)) {
    x[x %in% missingCodes] <- NA
    return(x)
}

add_paleon_grid <- function(data) {
    library(dplyr)
    points <- data %>% dplyr::select(x,y)
    coordinates(points) <- ~x+y
    proj4string(points) <- CRS('+init=epsg:3175')
    data <- data %>% mutate(cell = raster::extract(base_raster, points))
    return(data)
}

get_angle_inil <- function(bearings, degrees) {
    ##  This function converts the bearing and degree information
    ##  into 360-degree azimuth values.

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
        warning("get_angle_inil: Found non-zero degrees for cardinal bearings.")
        
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
  #  This function converts the azimuth strings to numeric, 360 degree values.
   angle <- matrix(as.numeric(NA), nrow(azimuth), ncol(azimuth))
  
  #  This is a special case, generally where the tree is plot center.
    angle[azimuth == '0'] <- 0

    # these will be assigned angle of 45
    azimuth <- gsub("XX", "  ", azimuth)
    azimuth <- gsub("--", "  ", azimuth)
    
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
        mutate(state = ifelse(grepl("Detroit", state), "MI", state)) %>%   
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
    ## don't treat distance of zero as being in same quad
    dist0 <- (!is.na(data$dist1) & data$dist1 == 0) | (!is.na(data$dist2) & data$dist2 == 0)
    same_quad <- same_quad & !dist0
    
    density <- rep(as.numeric(NA), nrow(data))

    ## fix upper radius searched and take 1 tree per that area as density
    max_dist_meters <- max_distance_surveyed * meters_per_link
    density[data$num_trees == 1] <- (1 / (pi * max_dist_meters^2)) * meters_sq_per_ha
    
    density[data$num_trees == 0] <- 0

    to_calc <- data$num_trees == 2
    sub <- data[to_calc, ]

    diam <- sub$diam1
    diam[is.na(diam)] <- 10
    ## distance in meters: distance plus one-half diameter of tree
    distm1 <- sub$dist1 * meters_per_link + 0.5 * diam * cm_per_inch / cm_per_m
    diam <- sub$diam2
    diam[is.na(diam)] <- 10
    ## distance in meters: distance plus one-half diameter of tree
    distm2 <- sub$dist2 * meters_per_link + 0.5 * diam * cm_per_inch / cm_per_m

    ##  From the formula,
    ##  lambda = kappa * theta * (q - 1)/(pi * n) * (q / sum_(1:q)(r^2))
    ##  here, n is equal to 1.
    ##  units are in stems / m^2
    q <- 2

    morisita <- ((q - 1) / pi) * (q / (distm1^2 + distm2^2)) * sub$full_corr 
    morisita <- morisita * meters_sq_per_ha

    density[to_calc] <- morisita

    cat("Found ", sum(same_quad), " corners with trees in same quadrant.\n")
    density[same_quad] <- NA

    density[density > max_density] <- max_density

    return(density)
}

calc_biomass_taxon <- function(num_trees, biomass1, biomass2, density, L3_tree1, L3_tree2, taxon) {
    biomass <- rep(0, length(num_trees))
    
    cond <- num_trees == 1 & L3_tree1 == taxon
    biomass[cond] <- biomass1[cond] * density[cond] / kg_per_Mg
    ## assign half biomass to taxon if either tree is the taxon
    cond <- num_trees == 2 & L3_tree1 == taxon 
    biomass[cond] <- biomass1[cond] * density[cond]  / (kg_per_Mg * 2)  ## 2 to account for taxon represents half the density
    ## if two trees of same taxon, the addition should handle this
    cond <- num_trees == 2 & L3_tree2 == taxon 
    biomass[cond] <- biomass[cond] +
        biomass2[cond] * density[cond]  / (kg_per_Mg * 2)

    ## handle case of two trees same taxon but one biomass is missing; use single biomass as the per-tree estimate
    cond <- num_trees == 2 & L3_tree1 == taxon & L3_tree2 == taxon & is.na(biomass1) & !is.na(biomass2)
    biomass[cond] <- biomass2[cond] * density[cond]  / kg_per_Mg
    cond <- num_trees == 2 & L3_tree1 == taxon & L3_tree2 == taxon & is.na(biomass2) & !is.na(biomass1)
    biomass[cond] <- biomass1[cond] * density[cond]  / kg_per_Mg
        
    return(biomass)
}

calc_basalarea_taxon <- function(num_trees, basalarea1, basalarea2, density, L3_tree1, L3_tree2, taxon) {
    basalarea <- rep(0, length(num_trees))
    
    cond <- num_trees == 1 & L3_tree1 == taxon
    basalarea[cond] <- basalarea1[cond] * density[cond] 
    ## assign half basalarea to taxon if either tree is the taxon
    cond <- num_trees == 2 & L3_tree1 == taxon 
    basalarea[cond] <- basalarea1[cond] * density[cond]  / 2  ## 2 to account for taxon represents half the density
    ## if two trees of same taxon, the addition should handle this
    cond <- num_trees == 2 & L3_tree2 == taxon 
    basalarea[cond] <- basalarea[cond] +
        basalarea2[cond] * density[cond]  / 2

    ## handle case of two trees same taxon but one basal area is missing; use single basal area as the per-tree estimate
    cond <- num_trees == 2 & L3_tree1 == taxon & L3_tree2 == taxon & is.na(basalarea1) & !is.na(basalarea2)
    basalarea[cond] <- basalarea2[cond] * density[cond]  
    cond <- num_trees == 2 & L3_tree1 == taxon & L3_tree2 == taxon & is.na(basalarea2) & !is.na(basalarea1)
    basalarea[cond] <- basalarea1[cond] * density[cond]  
        
    return(basalarea)
}

calc_density_taxon <- function(num_trees, density, L3_tree1, L3_tree2, taxon) {
    taxon_density <- rep(0, length(num_trees))
    
    cond <- num_trees == 1 & L3_tree1 == taxon
    taxon_density[cond] <- density[cond] 
    ## assign half biomass to taxon if either tree is the taxon
    cond <- num_trees == 2 & L3_tree1 == taxon 
    taxon_density[cond] <- density[cond]  / 2  ## 2 to account for taxon represents half the density
    ## if two trees of same taxon, the addition should handle this
    cond <- num_trees == 2 & L3_tree2 == taxon 
    taxon_density[cond] <- taxon_density[cond] + density[cond] / 2

    return(taxon_density)
}

wgt_mse <- function(n, y, yhat) {
    sum(n * (y - yhat)^2, na.rm = TRUE) / sum(n > 0, na.rm = TRUE)
}

calc_point_criterion <- function(pred_occ, pred_pot, n, y, mx, obj_fun = wgt_mse) {
    crit <- matrix(0, ncol(pred_occ), ncol(pred_pot))
    y[y > mx] <- mx
    for(i in seq_len(nrow(crit))) {
        for(j in seq_len(ncol(crit))) {
            tmp <- pred_occ[ , i] * pred_pot[ , j]
            tmp[tmp > mx] <- mx
            crit[i, j] <- obj_fun(n, y, tmp)
        }}
    dimnames(crit)[[1]] <- dimnames(pred_occ)[[2]]
    dimnames(crit)[[2]] <- dimnames(pred_pot)[[2]]
    return(crit)
}

calc_cov_criterion <- function(draws_logocc, draws_logpot, sig2, data, min_points = 60, n_draw = 250, seed = 1, type_pot = 'arith', scale = 1, size = 0.90) {
    set.seed(seed)
    if(is.null(min_points)) min_points <- 0
    wh <- which(data$points_total >= min_points)
    
    cov <- length <- loglength <- matrix(0, dim(draws_logocc)[[2]], dim(draws_logpot)[[2]])
    N <- nrow(data)

    if(scale) sig2 <- sig2 / (data$points_occ/scale)
    sig2[sig2 == Inf] <- 1e6  # just so calcs don't fail; these should be left out when [wh] applied
    for(i in seq_len(nrow(cov))) {
        for(j in seq_len(ncol(cov))) {
            tmp <- matrix(0, N, n_draw)
            for(k in seq_len(n_draw)) {
                yocc <- rbinom(N, data$points_total, exp(draws_logocc[ , i, k]))
                ypot <- rnorm(N, draws_logpot[ , j, k], sqrt(sig2[ , j]))
                if(type_pot == 'log_arith') ypot <- exp(ypot)
                tmp[ , k] <- ypot*yocc/data$points_total
            }
            qq <- apply(tmp, 1, quantile, c((1-size)/2, 1-(1-size)/2), na.rm = TRUE)
            qq <- qq[, wh]
            cov[i,j] <- mean(data$y[wh] < qq[2, ] & data$y[wh] > qq[1, ])
            length[i, j] <- median(qq[2, ] - qq[1, ])
            loglength[i, j] <- median(log(qq[2, ]) - log(qq[1, ]))
        }}
    dimnames(cov)[[1]] <- dimnames(length)[[1]] <- dimnames(loglength)[[1]] <- dimnames(draws_logocc)[[2]]
    dimnames(cov)[[2]] <- dimnames(length)[[2]] <- dimnames(loglength)[[2]] <- dimnames(draws_logpot)[[2]]
    return(list(cov = cov, length = length, loglength = loglength))

}

convert_chains_to_links <- function(dists) {
    ## x.5 decimals may well be links
    dec_valued <- round(dists) != dists
    dists_in_chains <- dec_valued &
        ( dists <= 1 | (dists > 1 & dists %% 1 != 0.5) )
    dists[dists_in_chains] <- dists[dists_in_chains] * links_per_chain
    return(dists)
}
