## Code for fitting occupancy and potential biomass/density models.

library(mgcv)

## CLEANUP: gamma doesn't respond to bam with fREML default

fit <- function(data, newdata, k_occ = NULL, k_pot = NULL, unc = FALSE, points_total = 'points_total',
                points_occ = 'points_occ', weight = points_occ, weight_scale = 1, avg = 'avg', geom_avg = 'geom_avg',
                gamma = 1, units = 'm', use_bam = FALSE, type_pot = 'arith', return_model = FALSE,
                save_draws = FALSE, num_draws = 250, bound_draws_low = FALSE, bound_draws_high = TRUE) {
    ## fits occ and pot components for single k values, with uncertainty if desired OR
    ## fits one or both components for one or more k values, without uncertainty
    
    if(!type_pot %in% c('arith', 'log_arith', 'geom'))
        stop("type_pot must be one of 'arith', 'log_arith', 'geom'")

    ## bam() handles big datasets more quickly and seems to be more numerically stable.
    if(use_bam)
        fitter <- bam else fitter <- gam
    
    if(units == 'm') {
        scaling <- 8000
    } else {
        if(units == 'km') {
            scaling <- 8
        } else stop("units should either be 'm' or 'km'")
    }
    
    newdata[ , c('x','y')] <- newdata[ , c('x','y')] / scaling
    data[ , c('x','y')] <- data[ , c('x','y')] / scaling 

    pred_occ <- pred_occ_se <- pred_pot <- pred <-
        draws <- draws_logocc <- draws_logocc_orig <- draws_logpot <- 
        model_occ <- model_pot <- NULL
    
  ###########################
  #  stage 1: fit occupancy #
  ###########################

    if(!is.null(k_occ)) {
        se.fit <- ifelse(length(k_occ) == 1, TRUE, FALSE)
        data$z <- cbind(data[[points_occ]], data[[points_total]] - data[[points_occ]])

        ## default bam() method is "fREML" which doesn't seem to respond to 'gamma'
        ## and is not the same method as used in gam()
        model_occ <- pred_occ <- list()
        for(k_idx in seq_along(k_occ)) {
            model_occ[[k_idx]] <- fitter(z ~ s(x, y, k = k_occ[k_idx]), data = data, family = 'binomial',
                                         gamma = gamma)
            pred_occ[[k_idx]] <- predict(model_occ[[k_idx]], newdata = newdata, type = 'response',
                                         se.fit = se.fit)
        }
    } else pred_occ <- 1
 
  ###################################
  #  stage 2: fit potential         #
  ###################################

    if(!is.null(k_pot)) {
        data <- data[data[[points_occ]] > 0, ]
        if(type_pot == 'arith') data$z <- data[[avg]]
        if(type_pot == 'log_arith') data$z <- log(data[[avg]])
        if(type_pot == 'geom') data$z <- data[[geom_avg]]
        data$weight <- data[[weight]] / weight_scale
        
        model_pot <- pred_pot <- list()
        for(k_idx in seq_along(k_pot)) {
            ## make sure k value not too large relative to sample size
            kval <- min(k_pot[k_idx], round(nrow(data)*0.9))

            model_pot[[k_idx]] <- fitter(z ~ s(x,y, k = kval), data = data, weights = weight,
                                     gamma = gamma)            
            pred_pot[[k_idx]] <- predict(model_pot[[k_idx]], newdata= newdata, type='response')
            if(type_pot != 'arith') {
                pred_pot[[k_idx]] <- exp(pred_pot[[k_idx]])
            } else pred_pot[[k_idx]][pred_pot[[k_idx]] < 0] <- 0
        }
    }
    
  #####################
  # predicted result  #
  #####################

    if((is.null(k_occ) || length(k_occ) == 1) && length(k_pot) == 1) {
        pred <- pred_occ * pred_pot
    } else pred <- 0   # don't do all pairs of possible preds across different k_occ and k_pot
  
    if(unc) {  # implement approx. Bayesian posterior draws following Wood (2004) section 4.8
        
        rmvn <- function(n, mu, sig) { ## MVN random deviates
            L <- mroot(sig); m <- ncol(L);
            t(mu + L %*% matrix(rnorm(m*n), m, n)) 
        }
        
        ## posterior draws of (log) occupancy
        if(!is.null(k_occ)) {
            draws_logocc <- array(0, c(nrow(newdata), length(k_occ), num_draws))
            dimnames(draws_logocc)[[2]] <- k_occ
            if(bound_draws_low || bound_draws_high) 
                draws_logocc_orig <- draws_logocc 
            for(k_idx in seq_along(k_occ)) {
                Xp <- predict(model_occ[[k_idx]], newdata = newdata, type="lpmatrix")
                draws_coef <- rmvn(num_draws , coef(model_occ[[k_idx]]), model_occ[[k_idx]]$Vp) 
                draws_linpred <- Xp %*% t(draws_coef)
                draws_logocc_tmp <- -log(1 + exp(-draws_linpred)) # log scale to add to log pot result
                if(bound_draws_low || bound_draws_high) {
                    draws_logocc_orig_tmp <- draws_logocc_tmp
                    ## Draws_linpred can have high variance near boundary,
                    ## where value of draws_linpred is very negative (so occ=0)
                    ## individual draws then can have high occ in those areas.
                    ## 
                    ## Draws_linpred can have high variance in small areas producing very positive
                    ## linpred values corresponding to very small pred_occ values.
                    ## 
                    ## Hacky fix that seems reasonable: draws bigger than 5x point estimate of occupancy
                    ## replaced with point prediction.
                    if(bound_draws_low) {
                        log_predocc_plus5 <- log(5) + log(pred_occ[[k_idx]])
                        draws_logocc_tmp <- apply(draws_logocc_tmp, 2, function(x) {
                            x[x > log_predocc_plus5] <- log(pred_occ[[k_idx]][x > log_predocc_plus5])
                            return(x)})
                    }
                    ## address numerical issue that seems to arise
                    ## (e.g., central MI in total FIA biomass)
                    ## that produces some draws where Pr(occ)=0 but point prediction is 1
                    if(bound_draws_high) 
                        draws_logocc_tmp[pred_occ[[k_idx]]> 0.999] <- 0
                } 
                draws_logocc[ , k_idx, ] <- draws_logocc_tmp
                if(bound_draws_low || bound_draws_high) 
                    draws_logocc_orig[ , k_idx, ] <- draws_logocc_orig_tmp
            }
        } 
        ## best to construct CIs on log scale and exponentiate endpoints
            
        ## posterior draws of (log) potential result
        if(!is.null(k_pot)) {
            draws_logpot <- array(0, c(nrow(newdata), length(k_pot), num_draws))
            dimnames(draws_logpot)[[2]] <- k_pot
            for(k_idx in seq_along(k_pot)) {
                Xp <- predict(model_pot[[k_idx]], newdata = newdata, type="lpmatrix")
                draws_coef <- rmvn(num_draws , coef(model_pot[[k_idx]]), model_pot[[k_idx]]$Vp) 
                draws_logpot_tmp <- Xp %*% t(draws_coef)
                if(type_pot == 'arith') {
                    draws_logpot_tmp[draws_logpot_tmp < 0] <- 0.01
                    draws_logpot_tmp <- log(draws_logpot_tmp)
                }
                draws_logpot[ , k_idx, ] <- draws_logpot_tmp
            }
        }
        
        if((is.null(k_occ) || length(k_occ) == 1) && length(k_pot) == 1) {
            if(is.null(k_occ)) {
                draws <- exp(draws_logpot[ , 1, ])
            } else draws <- exp(draws_logocc[ , 1, ] + draws_logpot[ , 1, ])
            pp.sd <- apply(draws, 1, sd)
            pred <- data.frame(mean = pred, sd = pp.sd)
        } 
        if(!save_draws)
            draws <- draws_logocc <- draws_logocc_orig <- draws_logpot <- NULL
    } 

    ## also reduce dim of draws_logocc, draws_logpot
    if(length(k_occ) == 1) {
        model_occ <- model_occ[[1]]
        tmpfit <- pred_occ[[1]]$fit
        pred_occ_se <- pred_occ[[1]]$se.fit
        pred_occ <- tmpfit
        if(unc && save_draws) {
            draws_logocc <- draws_logcc[[1]]
            if(bound_draws_low || bound_draws_high) 
                draws_logocc_orig <- draws_logocc_orig[[1]]
        }
    } else {
        if(length(k_occ) > 1) {
            names(model_occ) <- names(pred_occ) <- k_occ
            pred_occ <- as.matrix(as.data.frame(pred_occ))
            dimnames(pred_occ) <- NULL
            dimnames(pred_occ)[[2]] <- k_occ
        }
    }
    if(length(k_pot) == 1) {
        model_pot <- model_pot[[1]]
        pred_pot <- pred_pot[[1]]
        if(unc && save_draws)
            draws_logpot <- draws_logpot[[1]]
    } else {
        if(length(k_pot) > 1) {
            names(model_pot) <- names(pred_pot) <- k_pot
            pred_pot <- as.matrix(as.data.frame(pred_pot))
            dimnames(pred_pot) <- NULL
            dimnames(pred_pot)[[2]] <- k_pot
        }
    }

    if(!return_model) {
        model_occ <- NULL
        model_pot <- sapply(model_pot, '[[', 'sig2')
    }
    return(list(locs = data.frame(x = newdata$x*scaling, y = newdata$y*scaling),
                pred = pred, pred_occ = pred_occ, pred_occ_se = pred_occ_se,
                pred_pot = pred_pot, draws = draws,
                draws_logpot = draws_logpot, draws_logocc = draws_logocc,
                draws_logocc_orig = draws_logocc_orig,
                model_occ = model_occ, model_pot = model_pot,
                k_occ = k_occ, k_pot = k_pot))
}



