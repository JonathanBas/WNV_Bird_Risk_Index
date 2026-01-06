predic_mod <- function(resmod, nsim=1000) {
  post_expec_lin_pred2 <- resmod$summary.fitted.values[,"0.5quant"]

  post_expec_lin_pred <- as.vector(unlist(mclapply(X = resmod$marginals.fitted.values,
                                                   FUN = function(marg){
                                                     try_emarg <- try(inla.emarginal(fun = function(x){x},
                                                                                     marginal = marg), silent=T)
                                                     if(! is.numeric(try_emarg)){
                                                       try_emarg <- NA
                                                     }
                                                     try_emarg
                                                   })))
  post_expec_lin_pred[is.na(post_expec_lin_pred)] <- post_expec_lin_pred2[is.na(post_expec_lin_pred)]

  if(resmod$.args$family == "binomial"){
    vpred <- resmod$.args$Ntrials * post_expec_lin_pred
  }else if(resmod$.args$family == "poisson"){
    vpred <- post_expec_lin_pred
  }else if(resmod$.args$family == "nbinomial"){
    vpred <- post_expec_lin_pred
  }else if(resmod$.args$family == "zeroinflatedpoisson1"){
    vpred <- sapply(X = 1:length(post_expec_lin_pred),
                    FUN = function(id_i){
                      mean(rzip(n = nsim,
                                lambda = post_expec_lin_pred[id_i],
                                psi = resmod$summary.hyperpar["zero-probability parameter for zero-inflated poisson_1", "mean"]))})
  }else if(resmod$.args$family == "zeroinflatednbinomial1"){
    vpred <- sapply(X = 1:length(post_expec_lin_pred),
                    FUN = function(id_i){
                      mean(rzinb(n = nsim,
                                 mu = post_expec_lin_pred[id_i],
                                 size = resmod$summary.hyperpar["size for nbinomial zero-inflated observations", "mean"],
                                 psi = resmod$summary.hyperpar["zero-probability parameter for zero-inflated nbinomial_1", "mean"]))})
  }
  
  return(vpred)
}

comp_pit <- function(vec_obs, modres=NULL, vec_pred_poiss=NULL, which_id=NULL, n_breaks=10, tit_fig="Histogram of PIT", x_lab="Probability integral transform", y_lab="Density"){
  ## Adapted from "tscount" package
  ## Either "modres" (an INLA object) or directly the vector of predictions (from Poisson regression only) "vec_pred_poiss" is entered
  ## Either all the data (default), or possible to determine which IDs are used with "which_id" (example: training and testing datasets)
  
  breaks_hist <- seq(0, 1, length = n_breaks+1)
  hist_pit <- rep(0, length(breaks_hist))
  
  if(is.null(which_id)){
    n_obs_incl <- length(vec_obs)
    v_incl_obs <- 1:n_obs_incl
  }else{
    n_obs_incl <- length(which_id)
    v_incl_obs <- which_id
  }
  
  if(! is.null(modres)){
    post_expec_lin_pred2 <- modres$summary.fitted.values[,"0.5quant"]

    post_expec_lin_pred <- as.vector(unlist(mclapply(X = modres$marginals.fitted.values,
                                                     FUN = function(marg){
                                                       try_emarg <- try(inla.emarginal(fun = function(x){x},
                                                                                   marginal = marg), silent=T)
                                                       if(! is.numeric(try_emarg)){
                                                         try_emarg <- NA
                                                       }
                                                       try_emarg
                                                     })))
    post_expec_lin_pred[is.na(post_expec_lin_pred)] <- post_expec_lin_pred2[is.na(post_expec_lin_pred)]
  }
  
  p_func <- function(obser, id_k){
    if(is.null(modres)){
      ppois(q = obser, lambda = vec_pred_poiss[id_k])
    }else{
      
      if(modres$.args$family == "binomial"){
        pbinom(q = obser, size = modres$.args$Ntrials, prob = post_expec_lin_pred[id_k])
      }else if(modres$.args$family == "poisson"){
        ppois(q = obser, lambda = post_expec_lin_pred[id_k])
      }else if(modres$.args$family == "nbinomial"){
        pnbinom(q = obser,
                mu = post_expec_lin_pred[id_k],
                size = modres$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", "mean"])
      }else if(modres$.args$family == "zeroinflatedpoisson1"){
        pzip(q = obser,
             lambda = post_expec_lin_pred[id_k],
             psi = modres$summary.hyperpar["zero-probability parameter for zero-inflated poisson_1", "mean"])
      }else if(modres$.args$family == "zeroinflatednbinomial1"){
        pzinb(q = obser,
              mu = post_expec_lin_pred[id_k],
              size = modres$summary.hyperpar["size for nbinomial zero-inflated observations", "mean"],
              psi = modres$summary.hyperpar["zero-probability parameter for zero-inflated nbinomial_1", "mean"])
      }
    }
  }
  
  for(t in v_incl_obs){
    P_x <- p_func(vec_obs[t], t)
    
    if(vec_obs[t] != 0){
      P_x_1 <- p_func(vec_obs[t] - 1, t)
    }else{
      P_x_1 <- 0
    }
    hist_pit <- hist_pit + punif(q = breaks_hist, min = P_x_1, max = P_x) /n_obs_incl
  }
  
  histo <- list(breaks=breaks_hist, counts=diff(hist_pit)*n_obs_incl, density=diff(hist_pit), mids=(breaks_hist[-(n_breaks+1)]+breaks_hist[-1])/2, equidits=TRUE)
  class(histo) <- "histogram"
  plot(histo, main=tit_fig, xlab=x_lab, ylab=y_lab, freq=F)
  abline(h=1/n_breaks, lty="dashed", col="blue")
}

R2_mod <- function(res_mod, res_null_mod, type_R2="deviance", vec_subset=TRUE, dat_vect=NULL){
  if(is.null(dat_vect)){
    name_outcome <- as.character(res_mod$.args$formula)[2]
    v_data_outcome <- res_mod$.args$data[vec_subset, name_outcome]
  }else{
    v_data_outcome <- dat_vect[vec_subset]
  }
  
  v_modpred <- predic_mod(resmod=res_mod)[vec_subset]
  v_modnullpred <- predic_mod(resmod=res_null_mod)[vec_subset]

  if((res_mod$.args$family == "binomial") & (res_mod$.args$family == res_null_mod$.args$family)){
    N_trials <- res_mod$.args$Ntrials
    
    ## Definition of R2 by [Cox and Snell 1989], cited by [Nagelkerke 1991]
    LL_null <- sum(dbinom(x=v_data_outcome, size=N_trials, prob=v_modnullpred/N_trials, log=T))
    LL_mod <- sum(dbinom(x=v_data_outcome, size=N_trials, prob=v_modpred/N_trials, log=T))
    R2_CS <- 1 - exp(-(2/length(v_modpred)) * (LL_mod - LL_null))

    ## Definition of R2 by [Cameron & Windmeijer 1996]
    dev_null <- 2 * sum(dbinom(x=v_data_outcome, size=N_trials, prob=v_data_outcome/N_trials, log=T) - dbinom(x=v_data_outcome, size=N_trials, prob=v_modnullpred/N_trials, log=T)) # Definition of deviance from [Feng et al 2020]
    dev_mod <- 2 * sum(dbinom(x=v_data_outcome, size=N_trials, prob=v_data_outcome/N_trials, log=T) - dbinom(x=v_data_outcome, size=N_trials, prob=v_modpred/N_trials, log=T)) # Definition of deviance from [Feng et al 2020]
    R2_DEV <- 1 - dev_mod/dev_null
    
  }else if((res_mod$.args$family == "nbinomial") & (res_mod$.args$family == res_null_mod$.args$family)){
    param_size_mod <- res_mod$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", "mean"]
    param_size_null <- res_null_mod$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", "mean"]
    
    ## Definition of R2 by [Cox and Snell 1989], cited by [Nagelkerke 1991]
    LL_null <- sum(dnbinom(x=v_data_outcome, size=param_size_null, mu=v_modnullpred, log=T))
    LL_mod <- sum(dnbinom(x=v_data_outcome, size=param_size_mod, mu=v_modpred, log=T))
    R2_CS <- 1 - exp(-(2/length(v_modpred)) * (LL_mod - LL_null))
    
    ## Definition of R2 by [Cameron & Windmeijer 1996]
    dev_null <- 2 * sum(dnbinom(x=v_data_outcome, size=param_size_null, mu=v_data_outcome, log=T) - dnbinom(x=v_data_outcome, size=param_size_null, mu=v_modnullpred, log=T)) # Definition of deviance from [Feng et al 2020]
    dev_mod <- 2 * sum(dnbinom(x=v_data_outcome, size=param_size_mod, mu=v_data_outcome, log=T) - dnbinom(x=v_data_outcome, size=param_size_mod, mu=v_modpred, log=T)) # Definition of deviance from [Feng et al 2020]
    R2_DEV <- 1 - dev_mod/dev_null
    
  }else if((res_mod$.args$family == "poisson") & (res_mod$.args$family == res_null_mod$.args$family)){
    
    ## Definition of R2 by [Cox and Snell 1989], cited by [Nagelkerke 1991]
    LL_null <- sum(dpois(x=v_data_outcome, lambda=v_modnullpred, log=T))
    LL_mod <- sum(dpois(x=v_data_outcome, lambda=v_modpred, log=T))
    R2_CS <- 1 - exp(-(2/length(v_modpred)) * (LL_mod - LL_null))
    
    ## Definition of R2 by [Cameron & Windmeijer 1996]
    dev_null <- 2 * sum(dpois(x=v_data_outcome, lambda=v_data_outcome, log=T) - dpois(x=v_data_outcome, lambda=v_modnullpred, log=T)) # Definition of deviance from [Feng et al 2020]
    dev_mod <- 2 * sum(dpois(x=v_data_outcome, lambda=v_data_outcome, log=T) - dpois(x=v_data_outcome, lambda=v_modpred, log=T)) # Definition of deviance from [Feng et al 2020]
    R2_DEV <- 1 - dev_mod/dev_null
  }else{
    print("Likelihood distribution family not supported.")
    R2_CS = R2_DEV = NA
  }
  
  if(type_R2 == "deviance"){
    return(R2_DEV)
  }else if(type_R2 == "CS"){
    return(R2_CS)
  }
}

