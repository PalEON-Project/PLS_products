library(mgcv)

fit <- function(data, newdata, k_occ = NULL, k_pot = NULL, unc = FALSE, points_total = 'points_total', points_occ = 'points_occ', weight = points_occ, avg = 'avg', geom_avg = 'geom_avg', gamma = 1, units = 'm', use_bam = FALSE, type_pot = 'arith', return_model = FALSE, save_draws = FALSE, num_draws = 250) {
    ## fits occ and pot components for single k values, with uncertainty if desired OR
    ## fits one or both components for one or more k values, without uncertainty
    
    if(!type_pot %in% c('arith', 'log_arith', 'geom'))
        stop("type_pot must be one of 'arith', 'log_arith', 'geom'")

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

    pred_occ <- pred_pot <- pred <- draws <- model_occ <- model_pot <- NULL
    
  ###########################
  #  stage 1: fit occupancy #
  ###########################

    if(!is.null(k_occ)) {
        data$z <- cbind(data[[points_occ]], data[[points_total]] - data[[points_occ]])

        ## default method is "fREML" which doesn't seem to respond to 'gamma'
        ## and is not the same method as used in gam()
        model_occ <- pred_occ <- list()
        for(k_idx in seq_along(k_occ)) {
            model_occ[[k_idx]] <- fitter(z ~ s(x, y, k = k_occ[k_idx]), data = data, family = 'binomial',
                                     gamma = gamma, method = "GCV.Cp")
            pred_occ[[k_idx]] <- predict(model_occ[[k_idx]], newdata = newdata, type = 'response')
            
        }
        if(length(k_occ) == 1) {
            model_occ <- model_occ[[1]]
            pred_occ <- pred_occ[[1]]
        } else {
            names(model_occ) <- names(pred_occ) <- k_occ
            pred_occ <- as.matrix(as.data.frame(pred_occ))
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
        data$weight <- data[[weight]]
        
        model_pot <- pred_pot <- list()
        for(k_idx in seq_along(k_pot)) {
            model_pot[[k_idx]] <- fitter(z ~ s(x,y, k = k_pot[k_idx]), data = data, weights = weight,
                                     gamma = gamma, method = "GCV.Cp")            
            pred_pot[[k_idx]] <- predict(model_pot[[k_idx]], newdata= newdata, type='response')
            if(type_pot != 'arith') {
                pred_pot[[k_idx]] <- exp(pred_pot[[k_idx]])
            } else pred_pot[[k_idx]][pred_pot[[k_idx]] < 0] <- 0
        }
        if(length(k_pot) == 1) {
            model_pot <- model_pot[[1]]
            pred_pot <- pred_pot[[1]]
        } else {
            names(model_pot) <- names(pred_pot) <- k_pot
            pred_pot <- as.matrix(as.data.frame(pred_pot))
        }
        if(!return_model) model_pot <- NULL
    }
    
  #####################
  # predicted result  #
  #####################

    if((is.null(k_occ) || length(k_occ) == 1) && length(k_pot) == 1) {
        pred <- pred_occ * pred_pot
  
        if(unc) {  # implement approx. Bayesian posterior draws following Wood (2004) section 4.8

            rmvn <- function(n, mu, sig) { ## MVN random deviates
                L <- mroot(sig); m <- ncol(L);
                t(mu + L %*% matrix(rnorm(m*n), m, n)) 
            }
            
            ## posterior draws of (log) occupancy
            if(!is.null(k_occ)) {
                Xp <- predict(model_occ, newdata = newdata, type="lpmatrix")
                draws_coef <- rmvn(num_draws , coef(model_occ), model_occ$Vp) 
                draws_linpred <- Xp %*% t(draws_coef)
                draws_logocc <- draws_linpred - log(1 + exp(draws_linpred)) # log scale to add to log pot result
            } else draws_logocc <- 0
            ## best to construct CIs on log scale and exponentiate endpoints
            
            ## posterior draws of (log) potential result
            Xp <- predict(model_pot, newdata = newdata, type="lpmatrix")
            draws_coef <- rmvn(num_draws , coef(model_pot), model_pot$Vp) 
            draws_logpot <- Xp %*% t(draws_coef)
            if(type_pot == 'arith') {
                draws_logpot[draws_logpot < 0] <- 0.01
                draws_logpot <- log(draws_logpot)
            }
            
            draws <- exp(draws_logocc + draws_logpot)

            pp.sd <- apply(draws, 1, sd)
            if(!save_draws)
                draws <- NULL
            if(FALSE) {  ## not clear we want to do this
                ## zero out uncertainty outside of range boundary where weird things are happening
                pp.sd[pp.sd/result > 2 & result < 1] <- 0
                
                cat("Note that the Bayesian uncertainty in areas well outside of the range of a taxon has unreasonably large uncertainty, likely due to numerical issues in estimating very small probabilities. Uncertainty in these locations has been set to zero artificially.\n")
            }
            pred <- data.frame(mean = pred, sd = pp.sd)
        }
    } else if(unc) warning("more than one 'k' value -- not computing uncertainty.")
    if(!return_model) {
        model_occ <- NULL
        model_pot <- NULL
    }
    return(list(locs = data.frame(x = newdata$x*scaling, y = newdata$y*scaling),
                pred = pred, pred_occ = pred_occ, pred_pot = pred_pot, draws = draws,
                model_occ = model_occ, model_pot = model_pot,
                k_occ = k_occ, k_pot = k_pot))
}

fit_cv <- function(cell_full, k_occ, k_pot, n_cores) {
    pred_occ <- matrix(0, nrow(cell_full), length(k_occ))
    pred_pot_arith <- pred_pot_larith <- matrix(0, nrow(cell_full), length(k_pot))
    
    dimnames(pred_occ)[[2]] <- k_occ
    dimnames(pred_pot_arith)[[2]] <- dimnames(pred_pot_larith)[[2]] <- k_pot
    
    n_folds <- max(cell_full$fold)
    if(n_cores > 1) {
        library(doParallel)
        registerDoParallel(cores = n_cores)
        output <- foreach(i = seq_len(n_folds)) %dopar% {
            train <- cell_full %>% filter(fold != i)
            test <- cell_full %>% filter(fold == i)
            
            po <- fit(train, newdata = test, k_occ = k_occ, unc = FALSE)
            ppa <- fit(train, newdata = test, k_pot = k_pot, type_pot = 'arith', unc = FALSE)
            ppl <- fit(train, newdata = test, k_pot = k_pot, type_pot = 'log_arith', unc = FALSE)
            list(po, ppa, ppl)
            cat("n_fold: ", i, " ", date(), "\n")
        }
        for(i in seq_len(n_folds)) {
            pred_occ[cell_full$fold == i, ] <- output[[i]][[1]]$pred_occ
            pred_pot_arith[cell_full$fold == i, ] <- output[[i]][[2]]$pred_pot
            pred_pot_larith[cell_full$fold == i, ] <- output[[i]][[3]]$pred_pot
        }
    } else {
        for(i in seq_len(n_folds)) {
            train <- cell_full %>% filter(fold != i)
            test <- cell_full %>% filter(fold == i)
            
            po <- fit(train, newdata = test, k_occ = k_occ)
            ppa <- fit(train, newdata = test, k_pot = k_pot, type_pot = 'arith')
            ppl <- fit(train, newdata = test, k_pot = k_pot, type_pot = 'log_arith')
            
            pred_occ[cell_full$fold == i, ] <- po$pred_occ
            pred_pot_arith[cell_full$fold == i, ] <- ppa$pred_pot
            pred_pot_larith[cell_full$fold == i, ] <- ppl$pred_pot
            cat("n_fold: ", i, " ", date(), "\n")    
        }
    }
    return(list(pred_occ = pred_occ, pred_pot_arith = pred_pot_arith, pred_pot_larith = pred_pot_larith))
}






if(FALSE) { ## save for the moment in case needed later
    
    fit_occ <- function(data, newdata, k, total = 'total', count = 'count', gamma = 1, units = 'm', use_bam = FALSE, return_model = FALSE) {
        
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
        
        data$z <- cbind(data[[count]], data[[total]] - data[[count]])

        ## default method is "fREML" which doesn't seem to respond to 'gamma'
        ## and is not the same method as used in gam()
        model <- pr_occ <- list()
        for(k_idx in seq_along(k)) {
            model[[k_idx]] <- fitter(z ~ s(x, y, k = k[k_idx]), data = data, family = 'binomial',
                                     gamma = gamma, method = "GCV.Cp")
            pr_occ[[k_idx]] <- predict(model[[k_idx]], newdata = newdata, type = 'response')
            
        }
        if(length(k) == 1) {
            model <- model[[1]]
            pr_occ <- pr_occ[[1]]
        } else {
            names(model) <- names(pr_occ) <- k
            pr_occ <- as.matrix(as.data.frame(pr_occ))
        }
        if(!return_model) model <- NULL
        return(list(pr_occ = pr_occ, model = model))
    }
    
    fit_pot <- function(data, newdata, k = 100, total = 'total', count = 'count', weight = count, avg = 'avg', geom_avg = 'geom_avg', gamma = 1, units = 'm', use_bam = FALSE, type = 'arith', return_model = FALSE) {

        if(!type %in% c('arith', 'log_arith', 'geom'))
            stop("type must be one of 'arith', 'log_arith', 'geom'")
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
        
        data <- data[data[[count]] > 0, ]
        if(type == 'arith') data$z <- data[[avg]]
        if(type == 'log_arith') data$z <- log(data[[avg]])
        if(type == 'geom') data$z <- data[[geom_avg]]
        data$weight <- data[[weight]]
        
        model <- pr_pot <- list()
        for(k_idx in seq_along(k)) {
            model[[k_idx]] <- fitter(z ~ s(x,y, k = k[k_idx]), data = data, weights = weight,
                                     gamma = gamma, method = "GCV.Cp")
            
            pr_pot[[k_idx]] <- predict(model[[k_idx]], newdata= newdata, type='response')
            if(type != 'arith') {
                pr_pot[[k_idx]] <- exp(pr_pot[[k_idx]])
            } else pr_pot[[k_idx]][pr_pot[[k_idx]] < 0] <- 0
        }
        if(length(k) == 1) {
            model <- model[[1]]
            pr_pot <- pr_pot[[1]]
        } else {
            names(model) <- names(pr_pot) <- k
            pr_pot <- as.matrix(as.data.frame(pr_pot))
        }
        if(!return_model) model <- NULL
        return(list(pr_pot = pr_pot, model = model))
    }

}
