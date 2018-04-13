#' MCMC Sampling Algorithm for the Function-on-Scalars Regression Model
#' with pre-fixed basis functions
#'
#' Runs the MCMC for the function-on-scalars regression model based on
#' a knorwn expansion. Here, we assume the factor regression has independent errors.
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param use_fpca logical; if TRUE, use FPCA, otherwise use orthonormal B-splines
#' @param dist_reg_coef prior distribution for the regression coefficient evolution errors;
#' must be one of
#' \itemize{
#' \item "NIG" (Normal-inverse-Gamma)
#' \item "HS" (horseshoe prior)
#' \item "DHS" (dynamic horseshoe prior)
#' \item "BL" (Bayesian lasso)
#' \item "SV" (Gaussian stochastic volatility)
#' }
#' @param dist_reg_error prior distribution for the regression-level errors; same options as \code{dist_reg_coef}
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "sigma_et" (observation error SD; possibly dynamic)
#' \item "Yhat" (fitted values)
#' }
#' @param add_p_shrink logical; when TRUE, include a predictor-specific shrinkage term
#' (i.e., the scale term will be half-Cauchy)
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{Tm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x T x m},  may be inefficient
#'
#' @examples
#' # Fixme
#'
#' @import dsp KFAS truncdist
#' @export
fosr_basis = function(Y, tau, X = NULL, use_fpca = FALSE,
                dist_reg_coef = "HS", dist_reg_error = "NIG",
                nsave = 1000, nburn = 1000, nskip = 10,
                mcmc_params = list("beta", "fk", "alpha"),
                add_p_shrink = FALSE,
                computeDIC = TRUE){
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)
  
  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y); d = ncol(tau)
  
  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  #----------------------------------------------------------------------------
  # Initialize the main terms
  
  if(any(is.na(Y))) stop("Missing data not implemented for basis methods")
  Yna = Y # for consistency later
  
  if(use_fpca){
    if(d > 1) stop("FPCA only implemented for d = 1")
    
    Fmat = fpca.face(Y, center = TRUE, 
                     argvals = as.numeric(tau01),
                     knots=  max(20, min(ceiling(floor(median(rowSums(!is.na(Y))))/4), 150)),
                     pve=0.99)$efunctions
  } else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  #} else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(0.1*median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  
  # Initialize the factors
  Beta = tcrossprod(Y, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value
  
  # Initialize the conditional expectation:
  Btheta = tcrossprod(Beta, Fmat)
  
  # Initialize the (time-dependent) observation error SD:
  sigma_e = sd(Y - Btheta, na.rm=TRUE); sigma_et = rep(sigma_e, T)
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors (including the intercept)
  p = ncol(X)
  
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_tk = matrix(0, nrow = T, ncol = K) # Residuals
  for(k in 1:K) {
    fit_k = lm(Beta[,k] ~ X - 1)
    alpha_pk[,k] = fit_k$coef
    gamma_tk[,k] = Beta[,k] - fit_k$fitted
  }
  
  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])
  
  # Variance term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  # Initialize the corresponding SD term(s):
  evolParams_int = initEvolParams(omega = gamma_tk, evol_error = dist_reg_error)
  
  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
  
  # Update the error SD for gamma:
  sigma_gamma_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T)
  #----------------------------------------------------------------------------
  if(p > 1){
    omega = alpha_pk[-1,] # Not the intercept
    
    # First, initialize the variable selection scaling terms:
    if(add_p_shrink){
      lambda_j = apply(omega, 1, sd)
      lambda_0 = median(lambda_j) # Global term
      px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
      rep_lambda_j = rep(lambda_j, each = K) # Recurring term
    } else rep_lambda_j = rep(1, (p-1)*K) # Recurring scale term
    
    evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
    
    # This is the scaling SD (Note: doesn't actually depend on T)
    lambda_omega_pk = evolParams$sigma_wt
    
    # SD for omega:
    sigma_omega_pk = lambda_omega_pk*rep_lambda_j
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et = array(NA, c(nsave, T))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)
  
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    #----------------------------------------------------------------------------
    # Step 1: Impute the data, Y:
    #----------------------------------------------------------------------------
    # Not implemented in this case
    #if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Not necessary here
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    Y_tilde = tcrossprod(Y, t(Fmat)); sigma_tilde = sigma_et
    
    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 + sigma_gamma_tk[,k]^2)
      chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
      ell_k = crossprod(X, y_tilde_k/sigma_tilde_k^2)
      alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
    }
    
    # And sample the errors gamma_tk:
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + matrix(1/sigma_gamma_tk^2))
    postMean = matrix((Y_tilde - X%*%alpha_pk)/rep(sigma_tilde^2, times = K))*postSD^2
    gamma_tk = matrix(rnorm(n = T*K, mean = postMean, sd = postSD), nrow = T)
    
    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    # Not necessary, but update the conditional mean here
    Btheta = tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 5: Sample the observation error variance
    #----------------------------------------------------------------------------
    # Or use uniform prior?
    sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Btheta)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, T)
    
    #----------------------------------------------------------------------------
    # Step 6: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]
    
    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    a1_mu = uni.slice(a1_mu, g = function(a){
      dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_mu = uni.slice(a2_mu,g = function(a){
      sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(gamma_tk/evolParams_int$sigma_wt, ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
    a1_gamma = uni.slice(a1_gamma, g = function(a){
      dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_gamma = uni.slice(a2_gamma, g = function(a){
      sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Sample the corresponding prior variance term(s):
    evolParams_int = sampleEvolParams(omega = gamma_tk/rep(sigma_delta_k, each = T), evolParams =  evolParams_int, evol_error = dist_reg_error)
    
    # Update the error SD for gamma:
    sigma_gamma_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T)
    
    # Cap at machine epsilon:
    sigma_gamma_tk[which(sigma_gamma_tk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
    #----------------------------------------------------------------------------
    # Step 7: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      omega = alpha_pk[-1,] # Not the intercept
      
      #----------------------------------------------------------------------------
      # First part: j-specific shrinkage
      
      if(add_p_shrink){
        
        # Sum of squares (over k)
        ss_omega_j = rowSums((omega/lambda_omega_pk)^2)
        # Offset:
        ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
        
        # SD term:
        lambda_j = 1/sqrt(rgamma(n = p - 1,
                                 shape = 1/2 + K/2,
                                 rate = px_lambda_j + ss_omega_j/2))
        # Recurring term:
        rep_lambda_j = rep(lambda_j, each = K)
        
        # PX-j:
        px_lambda_j = rgamma(n = p - 1, shape = 1, rate = 1/lambda_j^2 + 1/lambda_0^2)
        
        # Global shrinkage:
        # SD term:
        lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_lambda_0 + sum(px_lambda_j)))
        
        # PX:
        px_lambda_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)
      }
      #----------------------------------------------------------------------------
      # Second part: jk-specific shrinkage
      
      # Sample the parameters
      evolParams = sampleEvolParams(omega = omega/rep_lambda_j, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
      
      # SD term:
      lambda_omega_pk = evolParams$sigma_wt
      #----------------------------------------------------------------------------
      # Last part:  update the variance parameters
      
      # SD for omega:
      sigma_omega_pk = lambda_omega_pk*rep_lambda_j
      
      # Cap at machine epsilon:
      sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    # Not necessary here
    #if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}
    
    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pk
        if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et[isave,] = sigma_et
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Btheta # + sigma_e*rnorm(length(Y))
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Btheta), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }
  
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('sigma_et', mcmc_params))) mcmc_output$sigma_et = post.sigma_et
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  
  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_et), m),
                            log = TRUE), na.rm=TRUE)
    
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d
    
    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))
  
  return (mcmc_output);
}
#' MCMC Sampling Algorithm for the Dynamic Function-on-Scalars Regression Model
#' with pre-fixed basis functions
#'
#' Runs the MCMC for the dynamic function-on-scalars regression model based on
#' a known expansion.
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param use_fpca logical; if TRUE, use FPCA, otherwise use orthonormal B-splines
#' @param factor_model model for the factor-level (regression) errors;
#' must be one of
#' \itemize{
#' \item "IND" (independent errors)
#' \item "AR" (stationary autoregression of order 1)
#' \item "RW" (random walk model)
#' }
#' @param use_dynamic_reg logical; if TRUE, regression coefficients are dynamic
#' (with random walk models), otherwise independent
#' @param dist_reg_coef prior distribution for the regression coefficient evolution errors;
#' must be one of
#' \itemize{
#' \item "NIG" (Normal-inverse-Gamma)
#' \item "HS" (horseshoe prior)
#' \item "DHS" (dynamic horseshoe prior)
#' \item "BL" (Bayesian lasso)
#' \item "SV" (Gaussian stochastic volatility)
#' }
#' @param dist_reg_error prior distribution for the regression-level errors; same options as \code{dist_reg_coef}
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept term for factor k)
#' \item "ar_phi" (AR coefficients for each k under AR(1) model)
#' \item "sigma_et" (observation error SD; possibly dynamic)
#' \item "Yhat" (fitted values)
#' }
#' @param add_p_shrink logical; when TRUE, include a predictor-specific shrinkage term
#' (i.e., the scale term will be half-Cauchy)
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model
#' for the observation error variance
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{Tm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x T x m},  may be inefficient
#'
#' @examples
#' # Fixme
#'
#' @import dsp KFAS truncdist
#' @export
dfosr_basis = function(Y, tau, X = NULL, use_fpca = FALSE,
                 factor_model = 'AR', use_dynamic_reg = TRUE,
                 dist_reg_coef = "DHS", dist_reg_error = "NIG",
                 nsave = 1000, nburn = 1000, nskip = 10,
                 mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                 add_p_shrink = FALSE,
                 use_obs_SV = FALSE,
                 computeDIC = TRUE){
  #----------------------------------------------------------------------------
  # Run some checks:
  # Convert to upper case, then check for matches to existing models:
  factor_model = toupper(factor_model); dist_reg_coef = toupper(dist_reg_coef); dist_reg_error = toupper(dist_reg_error)
  
  if(is.na(match(factor_model, c("RW", "AR", "IND"))))
    stop("The factor model must be one of 'RW', 'AR', or 'IND'")
  if(is.na(match(dist_reg_coef, c("DHS", "HS", "NIG", "BL"))))
    stop("The regression coefficient distribution must be one of 'DHS', 'HS', 'NIG', or 'BL'")
  if(is.na(match(dist_reg_error, c("DHS", "HS", "NIG", "BL"))))
    stop("The regression error distribution must be one of 'DHS', 'HS', 'NIG', or 'BL'")
  #----------------------------------------------------------------------------
  # Call the relevant function:
  
  if(factor_model == "RW")
    return(
      dfosr_rw_basis(Y = Y, tau = tau, X = X, use_fpca = use_fpca,
                     use_dynamic_reg = use_dynamic_reg,
                     dist_reg_coef = dist_reg_coef, dist_reg_error = dist_reg_error,
                     nsave = nsave, nburn = nburn, nskip = nskip,
                     mcmc_params = mcmc_params,
                     add_p_shrink = add_p_shrink,
                     use_obs_SV = use_obs_SV,
                     computeDIC = computeDIC)
    )
  
  if(factor_model == "AR")
    return(
      dfosr_ar_basis(Y = Y, tau = tau, X = X, use_fpca = use_fpca,
                     use_dynamic_reg = use_dynamic_reg,
                     dist_reg_coef = dist_reg_coef, dist_reg_error = dist_reg_error,
                     nsave = nsave, nburn = nburn, nskip = nskip,
                     mcmc_params = mcmc_params,
                     add_p_shrink = add_p_shrink,
                     use_obs_SV = use_obs_SV,
                     computeDIC = computeDIC)
    )
  
  if(factor_model == "IND")
    return(
      dfosr_ind_basis(Y = Y, tau = tau, X = X, use_fpca = use_fpca,
                      use_dynamic_reg = use_dynamic_reg,
                      dist_reg_coef = dist_reg_coef, dist_reg_error = dist_reg_error,
                      nsave = nsave, nburn = nburn, nskip = nskip,
                      mcmc_params = mcmc_params,
                      add_p_shrink = add_p_shrink,
                      use_obs_SV = use_obs_SV,
                      computeDIC = computeDIC)
    )
}
#' MCMC Sampling Algorithm for the Dynamic Function-on-Scalars Regression Model with Independent Errors
#'
#' Runs the MCMC for the dynamic function-on-scalars regression model based on
#' an FDLM-type expansion. Here, we assume the factor regression has independent errors.
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param use_fpca logical; if TRUE, use FPCA, otherwise use orthonormal B-splines
#' @param use_dynamic_reg logical; if TRUE, regression coefficients are dynamic
#' (with random walk models), otherwise independent
#' @param dist_reg_coef prior distribution for the regression coefficient evolution errors;
#' must be one of
#' \itemize{
#' \item "NIG" (Normal-inverse-Gamma)
#' \item "HS" (horseshoe prior)
#' \item "DHS" (dynamic horseshoe prior)
#' \item "BL" (Bayesian lasso)
#' \item "SV" (Gaussian stochastic volatility)
#' }
#' @param dist_reg_error prior distribution for the regression-level errors; same options as \code{dist_reg_coef}
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept term for factor k)
#' \item "sigma_et" (observation error SD; possibly dynamic)
#' \item "Yhat" (fitted values)
#' }
#' @param add_p_shrink logical; when TRUE, include a predictor-specific shrinkage term
#' (i.e., the scale term will be half-Cauchy)
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model
#' for the observation error variance
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{Tm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x T x m},  may be inefficient
#'
#' @examples
#' # Fixme
#'
#' @import dsp KFAS truncdist
#' @export
dfosr_ind_basis = function(Y, tau, X = NULL, use_fpca = FALSE,
                     use_dynamic_reg = TRUE,
                     dist_reg_coef = "DHS", dist_reg_error = "NIG",
                     nsave = 1000, nburn = 1000, nskip = 10,
                     mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                     add_p_shrink = FALSE,
                     use_obs_SV = FALSE,
                     computeDIC = TRUE){
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)
  
  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y); d = ncol(tau)
  
  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  #----------------------------------------------------------------------------
  # Initialize the main terms
  
  if(any(is.na(Y))) stop("Missing data not implemented for basis methods")
  Yna = Y # for consistency later
  
  if(use_fpca){
    if(d > 1) stop("FPCA only implemented for d = 1")
    
    Fmat = fpca.face(Y, center = TRUE, 
                     argvals = as.numeric(tau01),
                     knots=  max(20, min(ceiling(floor(median(rowSums(!is.na(Y))))/4), 150)),
                     pve=0.99)$efunctions
  } else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  #} else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(0.1*median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  
  # Initialize the factors
  Beta = tcrossprod(Y, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value

  # Initialize the conditional expectation:
  Btheta = tcrossprod(Beta, Fmat)
  
  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors:
  p = ncol(X)
  
  # This array is useful for setting up the DLM
  X.arr = array(0, c(K, K*p, T)); for(i in 1:T) X.arr[,,i] = blockDiag(matrix(X[i,], nrow=1) ,K)
  
  # Initialize the SSModel:
  kfas_model = update_kfas_model(Y.dlm = Beta, Zt = X.arr)
  
  # Indices identifying the intercept (or gamma_k term):
  ind.intercept = seq(1, K*p, by = p)
  
  # Identify the dynamic/non-dynamic components:
  if(!use_dynamic_reg) diag(kfas_model$R[,,1])[-ind.intercept] = 0
  #----------------------------------------------------------------------------
  # Overall mean term (and T x K case)
  mu_k = as.matrix(colMeans(Beta)); mu_tk = matrix(rep(mu_k, each =  T), nrow = T)
  
  # Variance term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  # Initialize the regression terms:
  
  # Evolution matrix (must be identity for static components)
  G_alpha = diag(K*p); diag(G_alpha)[ind.intercept] = 0
  
  # Update the SSModel object given the new parameters
  kfas_model = update_kfas_model(Y.dlm = Beta - mu_tk, Zt = X.arr,Gt = G_alpha, kfas_model = kfas_model)
  
  # Run the sampler
  alpha = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
  
  # Store as an array:
  alpha.arr = array(alpha, c(T, p, K))
  
  # Conditional mean from regression equation:
  for(k in 1:K) Beta[,k] = mu_k[k] + rowSums(X*alpha.arr[,,k])
  #----------------------------------------------------------------------------
  # Evolution error variance:
  Wt = array(diag(K*p), c(K*p, K*p, T)); W0 = diag(10^-4, K*p);
  
  # Intercept (or error) components:
  # use eta_tk for consistency w/ other functions
  eta_tk = as.matrix(alpha[,ind.intercept]) # gamma_tk = as.matrix(alpha[,ind.intercept])
  
  # Initialize the corresponding SD term(s):
  evolParams_int = initEvolParams(omega = eta_tk, evol_error = dist_reg_error)
  
  # MGP term:
  a1_eta = 2; a2_eta = 3;
  delta_eta_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
  
  # Update the error SD for gamma:
  sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T)
  
  # Evolution and initial variance:
  for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_eta_tk[-1, k]^2
  diag(W0[ind.intercept, ind.intercept]) = as.numeric(sigma_eta_tk[1,]^2)
  #----------------------------------------------------------------------------
  # Non-intercept term:
  if(p > 1){
    
    if(use_dynamic_reg){
      
      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])
      
      # Innovation:
      omega = diff(alpha_reg)
      
      # First, initialize the variable selection scaling terms:
      if(add_p_shrink){
        lambda_j = apply(array(omega, c(T-1, p-1, K)), 2, sd)
        lambda_0 = median(lambda_j) # Global term
        px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
        rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K) # Recurring term
      } else rep_lambda_j = rep(1, (T-1)*(p-1)*K) # Recurring scale term: set to 1
      
      # Initialize the variance components:
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))
      # This is the scaling SD
      lambda_omega_tpk = evolParams$sigma_wt
      
      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j
      
      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)
      
    } else{
      
      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      omega = t(alpha.arr[1, -1, ])
      
      # First, initialize the variable selection scaling terms:
      if(add_p_shrink){
        lambda_j = apply(omega, 2, sd)
        lambda_0 = median(lambda_j) # Global term
        px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
        rep_lambda_j = rep(lambda_j, each = K) # Recurring term
      } else rep_lambda_j = rep(1, (p-1)*K) # Recurring scale term
      
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      
      # This is the scaling SD (Note: doesn't actually depend on T)
      lambda_omega_pk = evolParams$sigma_wt
      
      # SD for omega:
      sigma_omega_pk = lambda_omega_pk*rep_lambda_j
      
      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))
    }
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, T, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et = array(NA, c(nsave, T))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)
  
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    #----------------------------------------------------------------------------
    # Step 1: Impute the data, Y:
    #----------------------------------------------------------------------------
    # Not implemented in this case
    #if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Not necessary here
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = tcrossprod(Y, t(Fmat)); sigma_tilde = sigma_et
    
    # Sanity check for Wt: if variances too large, KFAS will stop running
    Wt[which(Wt > 10^6, arr.ind = TRUE)] = 10^6; W0[which(W0 > 10^6, arr.ind = TRUE)] = 10^6
    
    # Update the SSModel object and sample:
    kfas_model = update_kfas_model(Y.dlm = Y_tilde - mu_tk,
                                   Zt = X.arr,
                                   sigma_et = sigma_tilde,
                                   Gt = G_alpha,
                                   Wt = Wt,
                                   W0 = W0,
                                   kfas_model = kfas_model)
    alpha = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
    
    # Store as an array:
    alpha.arr = array(alpha, c(T, p, K))
    
    # Update the factors:
    for(k in 1:K) Beta[,k] = mu_k[k] + rowSums(X*alpha.arr[,,k])
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    # Not necessary, but update the conditional mean here
    Btheta = tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 5: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Btheta, svParams)
      sigma_et = svParams$sigma_et
    } else {
      # Or use uniform prior?
      sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Btheta)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, T)
    }
    #----------------------------------------------------------------------------
    # Step 6: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    # Centerend and non-centered:
    gamma_tk =  as.matrix(alpha[,ind.intercept])
    gamma_tk_c = gamma_tk + mu_tk
    
    # Sample the unconditional mean term:
    postSD = 1/sqrt(colSums(1/sigma_eta_tk^2) + 1/sigma_mu_k^2)
    postMean = (colSums(gamma_tk_c/sigma_eta_tk^2))*postSD^2
    mu_k = rnorm(n = K, mean = postMean, sd = postSD)
    mu_tk = matrix(rep(mu_k, each =  T), nrow = T)
    
    # And update the non-centered parameter:
    gamma_tk = gamma_tk_c - mu_tk
    
    # This is the parameter we use (for consistency)
    eta_tk = gamma_tk
    
    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    a1_mu = uni.slice(a1_mu, g = function(a){
      dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_mu = uni.slice(a2_mu,g = function(a){
      sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                            delta.h = delta_eta_k,
                            a1 = a1_eta, a2 = a2_eta)
    sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
    a1_eta = uni.slice(a1_eta, g = function(a){
      dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_eta = uni.slice(a2_eta, g = function(a){
      sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Sample the corresponding prior variance term(s):
    evolParams_int = sampleEvolParams(omega = eta_tk/rep(sigma_delta_k, each = T), evolParams =  evolParams_int, evol_error = dist_reg_error)
    
    # Update the error SD for gamma:
    sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T)
    
    # Cap at machine epsilon:
    sigma_eta_tk[which(sigma_eta_tk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
    
    # Evolution and initial variance:
    for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_eta_tk[-1, k]^2
    diag(W0[ind.intercept, ind.intercept]) = as.numeric(sigma_eta_tk[1,]^2)
    #----------------------------------------------------------------------------
    # Step 7: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      
      if(use_dynamic_reg){
        
        # Dynamic setting
        
        # Regression (non-intercept) coefficients
        alpha_reg = as.matrix(alpha[,-ind.intercept])
        
        # Random walk, so compute difference for innovations:
        omega = diff(alpha_reg)
        
        #----------------------------------------------------------------------------
        if(add_p_shrink){
          # First part: j-specific shrinkage
          # Scale omega by lambda_omega_tpk:
          omega_scale = omega/lambda_omega_tpk
          # Sum of squares (over t,k)
          ss_omega_j = apply(array(omega_scale^2, c(T-1, p-1, K)), 2, sum)
          # Offset:
          ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
          
          # SD term:
          lambda_j = 1/sqrt(rgamma(n = p - 1,
                                   shape = 1/2 + (T-1)*K/2,
                                   rate = px_lambda_j + ss_omega_j/2))
          
          # Sample the PX, global parameters outside this if-statement
          
          # Recurring term:
          rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K)
        }
        #----------------------------------------------------------------------------
        # Second part: tjk-specific shrinkage
        
        # Scale by lambda_j
        omega_scale = omega/rep_lambda_j
        
        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)
        
        # SD term:
        lambda_omega_tpk = evolParams$sigma_wt
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters
        
        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j
        
        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
        
        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)
        
      } else{
        
        # Non-dynamic setting: grab the first one (all the same),
        # and store as (K x p-1) matrix
        omega = t(alpha.arr[1, -1, ])
        
        #----------------------------------------------------------------------------
        if(add_p_shrink){
          # First part: j-specific shrinkage
          
          # Sum of squares (over k)
          ss_omega_j = colSums((omega/lambda_omega_pk)^2)
          # Offset:
          ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
          
          # SD term:
          lambda_j = 1/sqrt(rgamma(n = p - 1,
                                   shape = 1/2 + K/2,
                                   rate = px_lambda_j + ss_omega_j/2))
          # Recurring term:
          rep_lambda_j = rep(lambda_j, each = K)
        }
        #----------------------------------------------------------------------------
        # Second part: jk-specific shrinkage
        
        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega/rep_lambda_j, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        
        # SD term:
        lambda_omega_pk = evolParams$sigma_wt
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters
        
        # SD for omega:
        sigma_omega_pk = lambda_omega_pk*rep_lambda_j
        
        # Cap at machine epsilon:
        sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
        
        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))
      }
      #----------------------------------------------------------------------------
      if(add_p_shrink){
        # In both dynamic/non-dynamic cases, need to sample the PX-parameters
        # for lambda_j and global shrinkage parameters
        
        # PX-j:
        px_lambda_j = rgamma(n = p - 1, shape = 1, rate = 1/lambda_j^2 + 1/lambda_0^2)
        
        # Global shrinkage:
        # SD term:
        lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_lambda_0 + sum(px_lambda_j)))
        
        # PX:
        px_lambda_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    # Not necessary here
    #if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}
    
    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,,] = alpha.arr
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et[isave,] = sigma_et
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Btheta # + sigma_e*rnorm(length(Y))
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Btheta), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }
  
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_et', mcmc_params))) mcmc_output$sigma_et = post.sigma_et
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  
  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_et), m),
                            log = TRUE), na.rm=TRUE)
    
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d
    
    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))
  
  return (mcmc_output);
}
#' MCMC Sampling Algorithm for the Dynamic Function-on-Scalars Regression Model with AR(1) Errors
#'
#' Runs the MCMC for the dynamic function-on-scalars regression model based on
#' an FDLM-type expansion. Here, we assume the factor regression has AR(1) errors.
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param use_fpca logical; if TRUE, use FPCA, otherwise use orthonormal B-splines
#' @param use_dynamic_reg logical; if TRUE, regression coefficients are dynamic
#' (with random walk models), otherwise independent
#' @param dist_reg_coef prior distribution for the regression coefficient evolution errors;
#' must be one of
#' \itemize{
#' \item "NIG" (Normal-inverse-Gamma)
#' \item "HS" (horseshoe prior)
#' \item "DHS" (dynamic horseshoe prior)
#' \item "BL" (Bayesian lasso)
#' \item "SV" (Gaussian stochastic volatility)
#' }
#' @param dist_reg_error prior distribution for the regression-level errors; same options as \code{dist_reg_coef}
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept term for factor k)
#' \item "ar_phi" (AR coefficients for each k under AR(1) model)
#' \item "sigma_et" (observation error SD; possibly dynamic)
#' \item "Yhat" (fitted values)
#' }
#' @param add_p_shrink logical; when TRUE, include a predictor-specific shrinkage term
#' (i.e., the scale term will be half-Cauchy)
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model
#' for the observation error variance
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{Tm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x T x m},  may be inefficient
#'
#' @examples
#' # Fixme
#'
#' @import dsp KFAS truncdist
#' @export
dfosr_ar_basis = function(Y, tau, X = NULL, use_fpca = FALSE,
                    use_dynamic_reg = TRUE,
                    dist_reg_coef = "DHS", dist_reg_error = "NIG",
                    nsave = 1000, nburn = 1000, nskip = 10,
                    mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                    add_p_shrink = FALSE,
                    use_obs_SV = FALSE,
                    computeDIC = TRUE){
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)
  
  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y); d = ncol(tau)
  
  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  
  if(any(is.na(Y))) stop("Missing data not implemented for basis methods")
  Yna = Y # for consistency later
  
  if(use_fpca){
    if(d > 1) stop("FPCA only implemented for d = 1")
    
    Fmat = fpca.face(Y, center = TRUE, 
                     argvals = as.numeric(tau01),
                     knots=  max(20, min(ceiling(floor(median(rowSums(!is.na(Y))))/4), 150)),
                     pve=0.99)$efunctions
  } else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  #} else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(0.1*median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  
  # Initialize the factors
  Beta = tcrossprod(Y, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value
  
  # Initialize the conditional expectation:
  Btheta = tcrossprod(Beta, Fmat)
  
  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors:
  p = ncol(X)
  
  # This array is useful for setting up the DLM
  X.arr = array(0, c(K, K*p, T)); for(i in 1:T) X.arr[,,i] = blockDiag(matrix(X[i,], nrow=1) ,K)
  
  # Initialize the SSModel:
  kfas_model = update_kfas_model(Y.dlm = Beta, Zt = X.arr)
  
  # Indices identifying the intercept (or gamma_k term):
  ind.intercept = seq(1, K*p, by = p)
  
  # Identify the dynamic/non-dynamic components:
  if(!use_dynamic_reg) diag(kfas_model$R[,,1])[-ind.intercept] = 0
  #----------------------------------------------------------------------------
  # Overall mean term (and T x K case)
  mu_k = as.matrix(colMeans(Beta)); mu_tk = matrix(rep(mu_k, each =  T), nrow = T)
  
  # Variance term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  # AR(1) Evolution Matrix
  
  # Evolution matrix (must be identity for static components)
  G_alpha = diag(K*p)
  
  # AR(1) coefficients:
  ar_int = apply(Beta - mu_tk, 2, function(x) lm(x[-1] ~ - 1 +  x[-length(x)])$coef)
  
  # Stationarity fix:
  ar_int[which(abs(ar_int) > 0.95)] = 0.8*sign(ar_int[which(abs(ar_int) > 0.95)])
  
  # Update the matrix:
  diag(G_alpha)[ind.intercept] = ar_int
  #----------------------------------------------------------------------------
  # Initialize the regression terms:
  
  # Update the SSModel object given the new parameters
  kfas_model = update_kfas_model(Y.dlm = Beta - mu_tk, Zt = X.arr,Gt = G_alpha, kfas_model = kfas_model)
  
  # Run the sampler
  alpha = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
  
  # Store as an array:
  alpha.arr = array(alpha, c(T, p, K))
  
  # Conditional mean from regression equation:
  for(k in 1:K) Beta[,k] = mu_k[k] + rowSums(X*alpha.arr[,,k])
  #----------------------------------------------------------------------------
  # Evolution error variance:
  Wt = array(diag(K*p), c(K*p, K*p, T)); W0 = diag(10^-4, K*p);
  
  # Intercept (or gamma) components:
  gamma_tk =  as.matrix(alpha[,ind.intercept])
  
  # Then subtract the AR(1) part:
  eta_tk = gamma_tk[-1,] -  t(ar_int*t(gamma_tk[-T,]))
  
  # Initialize the corresponding prior variance term(s):
  evolParams_int = initEvolParams(omega = eta_tk, evol_error = dist_reg_error)
  evolParams0_int = initEvol0(mu0 = as.matrix(gamma_tk[1,]))
  
  # MGP term:
  a1_eta = 2; a2_eta = 3;
  delta_eta_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
  
  # Update the error SD for gamma:
  sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T-1)
  
  # Evolution and initial variance:
  for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_eta_tk[,k]^2
  diag(W0[ind.intercept, ind.intercept]) = as.numeric(evolParams0_int$sigma_w0^2)
  #----------------------------------------------------------------------------
  # Non-intercept term:
  if(p > 1){
    
    if(use_dynamic_reg){
      
      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])
      
      # Innovation:
      omega = diff(alpha_reg)
      
      # First, initialize the variable selection scaling terms:
      if(add_p_shrink){
        lambda_j = apply(array(omega, c(T-1, p-1, K)), 2, sd)
        lambda_0 = median(lambda_j) # Global term
        px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
        rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K) # Recurring term
      } else rep_lambda_j = rep(1, (T-1)*(p-1)*K) # Recurring scale term: set to 1
      
      # Initialize the variance components:
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))
      # This is the scaling SD
      lambda_omega_tpk = evolParams$sigma_wt
      
      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j
      
      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)
      
    } else{
      
      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      omega = t(alpha.arr[1, -1, ])
      
      # First, initialize the variable selection scaling terms:
      if(add_p_shrink){
        lambda_j = apply(omega, 2, sd)
        lambda_0 = median(lambda_j) # Global term
        px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
        rep_lambda_j = rep(lambda_j, each = K) # Recurring term
      } else rep_lambda_j = rep(1, (p-1)*K) # Recurring scale term
      
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      
      # This is the scaling SD (Note: doesn't actually depend on T)
      lambda_omega_pk = evolParams$sigma_wt
      
      # SD for omega:
      sigma_omega_pk = lambda_omega_pk*rep_lambda_j
      
      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))
    }
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, T, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et = array(NA, c(nsave, T))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, T, m))
  if(!is.na(match('ar_phi', mcmc_params))) post.ar_phi = array(NA, c(nsave, K))
  if(computeDIC) post_loglike = numeric(nsave)
  
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    #----------------------------------------------------------------------------
    # Step 1: Impute the data, Y:
    #----------------------------------------------------------------------------
    # Not implemented in this case
    #if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Not necessary here
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = tcrossprod(Y, t(Fmat)); sigma_tilde = sigma_et
    
    # Sanity check for Wt: if variances too large, KFAS will stop running
    Wt[which(Wt > 10^6, arr.ind = TRUE)] = 10^6; W0[which(W0 > 10^6, arr.ind = TRUE)] = 10^6
    
    # Update the SSModel object and sample:
    kfas_model = update_kfas_model(Y.dlm = Y_tilde - mu_tk,
                                   Zt = X.arr,
                                   sigma_et = sigma_tilde,
                                   Gt = G_alpha,
                                   Wt = Wt,
                                   W0 = W0,
                                   kfas_model = kfas_model)
    alpha = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
    
    # Store as an array:
    alpha.arr = array(alpha, c(T, p, K))
    
    # Update the factors:
    for(k in 1:K) Beta[,k] = mu_k[k] + rowSums(X*alpha.arr[,,k])
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    # Not necessary, but update the conditional mean here
    Btheta = tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 5: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Btheta, svParams)
      sigma_et = svParams$sigma_et
    } else {
      # Or use uniform prior?
      sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Btheta)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, T)
    }
    #----------------------------------------------------------------------------
    # Step 6: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    # Centerend and non-centered:
    gamma_tk =  as.matrix(alpha[,ind.intercept])
    gamma_tk_c = gamma_tk + mu_tk
    
    # Sample the unconditional mean term:
    mu_k = sampleARmu(yt = gamma_tk_c,
                      phi_j = ar_int,
                      sigma_tj = sigma_eta_tk,
                      priorPrec = 1/sigma_mu_k^2)
    mu_tk = matrix(rep(mu_k, each =  T), nrow = T)
    
    # And update the non-centered parameter:
    gamma_tk = gamma_tk_c - mu_tk
    
    # AR(1) coefficients:
    ar_int = sampleARphi(yt = gamma_tk,
                         phi_j = ar_int,
                         sigma_tj = sigma_eta_tk,
                         prior_phi = c(5,2))
    #prior_phi = NULL)
    diag(G_alpha)[ind.intercept] = ar_int
    
    # Then subtract the AR(1) part:
    eta_tk = gamma_tk[-1,] -  t(ar_int*t(gamma_tk[-T,]))
    
    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    a1_mu = uni.slice(a1_mu, g = function(a){
      dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_mu = uni.slice(a2_mu,g = function(a){
      sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                            delta.h = delta_eta_k,
                            a1 = a1_eta, a2 = a2_eta)
    sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
    a1_eta = uni.slice(a1_eta, g = function(a){
      dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_eta = uni.slice(a2_eta, g = function(a){
      sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Sample the corresponding prior variance term(s):
    evolParams_int = sampleEvolParams(omega = eta_tk/rep(sigma_delta_k, each = T-1), evolParams =  evolParams_int, evol_error = dist_reg_error)
    evolParams0_int = sampleEvol0(mu0 = as.matrix(gamma_tk[1,]), evolParams0 = evolParams0_int)
    
    # Update the error SD for gamma:
    sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T-1)
    
    # Cap at machine epsilon:
    sigma_eta_tk[which(sigma_eta_tk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
    
    # Evolution and initial variance:
    for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_eta_tk[,k]^2
    diag(W0[ind.intercept, ind.intercept]) = as.numeric(evolParams0_int$sigma_w0^2)
    #----------------------------------------------------------------------------
    # Step 7: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      
      if(use_dynamic_reg){
        
        # Dynamic setting
        
        # Regression (non-intercept) coefficients
        alpha_reg = as.matrix(alpha[,-ind.intercept])
        
        # Random walk, so compute difference for innovations:
        omega = diff(alpha_reg)
        
        #----------------------------------------------------------------------------
        if(add_p_shrink){
          # First part: j-specific shrinkage
          # Scale omega by lambda_omega_tpk:
          omega_scale = omega/lambda_omega_tpk
          # Sum of squares (over t,k)
          ss_omega_j = apply(array(omega_scale^2, c(T-1, p-1, K)), 2, sum)
          # Offset:
          ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
          
          # SD term:
          lambda_j = 1/sqrt(rgamma(n = p - 1,
                                   shape = 1/2 + (T-1)*K/2,
                                   rate = px_lambda_j + ss_omega_j/2))
          
          # Sample the PX, global parameters outside this if-statement
          
          # Recurring term:
          rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K)
        }
        #----------------------------------------------------------------------------
        # Second part: tjk-specific shrinkage
        
        # Scale by lambda_j
        omega_scale = omega/rep_lambda_j
        
        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)
        
        # SD term:
        lambda_omega_tpk = evolParams$sigma_wt
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters
        
        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j
        
        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
        
        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)
        
      } else{
        
        # Non-dynamic setting: grab the first one (all the same),
        # and store as (K x p-1) matrix
        omega = t(alpha.arr[1, -1, ])
        
        #----------------------------------------------------------------------------
        if(add_p_shrink){
          # First part: j-specific shrinkage
          
          # Sum of squares (over k)
          ss_omega_j = colSums((omega/lambda_omega_pk)^2)
          # Offset:
          ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
          
          # SD term:
          lambda_j = 1/sqrt(rgamma(n = p - 1,
                                   shape = 1/2 + K/2,
                                   rate = px_lambda_j + ss_omega_j/2))
          # Recurring term:
          rep_lambda_j = rep(lambda_j, each = K)
        }
        #----------------------------------------------------------------------------
        # Second part: jk-specific shrinkage
        
        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega/rep_lambda_j, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        
        # SD term:
        lambda_omega_pk = evolParams$sigma_wt
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters
        
        # SD for omega:
        sigma_omega_pk = lambda_omega_pk*rep_lambda_j
        
        # Cap at machine epsilon:
        sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
        
        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))
      }
      #----------------------------------------------------------------------------
      if(add_p_shrink){
        # In both dynamic/non-dynamic cases, need to sample the PX-parameters
        # for lambda_j and global shrinkage parameters
        
        # PX-j:
        px_lambda_j = rgamma(n = p - 1, shape = 1, rate = 1/lambda_j^2 + 1/lambda_0^2)
        
        # Global shrinkage:
        # SD term:
        lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_lambda_0 + sum(px_lambda_j)))
        
        # PX:
        px_lambda_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    # Not necessary here
    #if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}
    
    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,,] = alpha.arr
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et[isave,] = sigma_et
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Btheta # + sigma_e*rnorm(length(Y))
        if(!is.na(match('ar_phi', mcmc_params))) post.ar_phi[isave,] = ar_int
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Btheta), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }
  
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_et', mcmc_params))) mcmc_output$sigma_et = post.sigma_et
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('ar_phi', mcmc_params))) mcmc_output$ar_phi = post.ar_phi
  
  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_et), m),
                            log = TRUE), na.rm=TRUE)
    
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d
    
    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))
  
  return (mcmc_output);
}
dfosr_basis_ar_k = function(Y, tau, X = NULL, use_fpca = FALSE,
                            use_dynamic_reg = TRUE,
                            dist_reg_coef = "NIG", dist_reg_error = "NIG",
                            nsave = 1000, nburn = 1000, nskip = 3,
                            mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                            add_p_shrink = FALSE,
                            use_obs_SV = FALSE,
                            computeDIC = TRUE){
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)
  
  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y); d = ncol(tau)
  
  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  
  if(any(is.na(Y))) stop("Missing data not implemented for basis methods")
  Yna = Y # for consistency later
  
  if(use_fpca){
    if(d > 1) stop("FPCA only implemented for d = 1")
    
    Fmat = fpca.face(Y, center = TRUE, 
                     argvals = as.numeric(tau01),
                     knots=  max(15, min(ceiling(floor(median(rowSums(!is.na(Y))))/4), 150)),
                     pve=0.99)$efunctions
  } else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  #} else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(0.1*median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  
  # Initialize the factors
  Beta = tcrossprod(Y, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value
  
  # Initialize the conditional expectation:
  Btheta = tcrossprod(Beta, Fmat)
  
  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors:
  p = ncol(X)
  
  # Initialize the SSModel:
  X.arr = array(t(X), c(1, p, T))
  kfas_model = update_kfas_model(Y.dlm = as.matrix(Beta[,1]), Zt = X.arr)
  
  # Identify the dynamic/non-dynamic components:
  if(p > 1 && !use_dynamic_reg) diag(kfas_model$R[,,1])[-1] = 0
  #----------------------------------------------------------------------------
  # Overall mean term (and T x K case)
  mu_k = as.matrix(colMeans(Beta)); mu_tk = matrix(rep(mu_k, each =  T), nrow = T)
  
  # Variance term for mu_k:
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  # AR(1) Evolution Matrix
  G_alpha = diag(p) # Replace the intercept terms as needed
  
  # AR(1) coefficients:
  ar_int = apply(Beta - mu_tk, 2, function(x) lm(x[-1] ~ - 1 +  x[-length(x)])$coef)
  
  # Stationarity fix:
  ar_int[which(abs(ar_int) > 0.95)] = 0.8*sign(ar_int[which(abs(ar_int) > 0.95)])
  
  #----------------------------------------------------------------------------
  # Initialize the regression terms:
  alpha.arr = array(0, c(T, p, K))
  for(k in 1:K){
    # Update the evoluation matrix
    G_alpha[1,1] = ar_int[k]
    
    # Update the SSModel object given the new parameters
    kfas_model = update_kfas_model(Y.dlm = as.matrix(Beta[,k] - mu_k[k]),
                                   Zt = X.arr,
                                   Gt = G_alpha, 
                                   kfas_model = kfas_model)
    # Run the sampler
    alpha.arr[,,k] = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
    
    # Conditional mean from regression equation:
    Beta[,k] = mu_k[k] + rowSums(X*alpha.arr[,,k])
  }
  
  #----------------------------------------------------------------------------
  # Evolution error variance:
  Wt = array(diag(p), c(p, p, T)); W0 = diag(10^-4, p);
  
  # Intercept (or gamma) components:
  gamma_tk =  matrix(alpha.arr[,1,], nrow = T)
  
  # Then subtract the AR(1) part:
  eta_tk = gamma_tk[-1,] -  t(ar_int*t(gamma_tk[-T,]))
  
  # Initialize the corresponding prior variance term(s):
  evolParams_int = initEvolParams(omega = eta_tk, evol_error = dist_reg_error)
  evolParams0_int = initEvol0(mu0 = as.matrix(gamma_tk[1,]))
  
  # MGP term:
  a1_eta = 2; a2_eta = 3;
  delta_eta_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
  
  # Update the error SD for gamma:
  sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T-1)
  #----------------------------------------------------------------------------
  # Non-intercept term:
  if(p > 1){
    
    if(use_dynamic_reg){
      
      # Dynamic setting:
      alpha_reg = matrix(alpha.arr[,-1,], nrow = T)
      
      # Innovation:
      omega = diff(alpha_reg)
      
      # First, initialize the variable selection scaling terms:
      if(add_p_shrink){
        lambda_j = apply(array(omega, c(T-1, p-1, K)), 2, sd)
        lambda_0 = median(lambda_j) # Global term
        px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
        rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K) # Recurring term
      } else rep_lambda_j = rep(1, (T-1)*(p-1)*K) # Recurring scale term: set to 1
      
      # Initialize the variance components:
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))
      # This is the scaling SD
      lambda_omega_tpk = evolParams$sigma_wt
      
      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j
      
    } else{
      
      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      #omega = t(alpha.arr[1, -1, ])
      omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)
      
      # First, initialize the variable selection scaling terms:
      if(add_p_shrink){
        lambda_j = apply(omega, 2, sd)
        lambda_0 = median(lambda_j) # Global term
        px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
        rep_lambda_j = rep(lambda_j, each = K) # Recurring term
      } else rep_lambda_j = rep(1, (p-1)*K) # Recurring scale term
      
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      
      # This is the scaling SD (Note: doesn't actually depend on T)
      lambda_omega_kp = evolParams$sigma_wt
      
      # SD for omega:
      sigma_omega_kp = lambda_omega_kp*rep_lambda_j
    }
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, T, p, K))
  if(!is.na(match('mu_k', mcmc_params))) post.mu_k = array(NA, c(nsave, K))
  if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et = array(NA, c(nsave, T))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, T, m))
  if(!is.na(match('ar_phi', mcmc_params))) post.ar_phi = array(NA, c(nsave, K))
  if(computeDIC) post_loglike = numeric(nsave)
  
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    #----------------------------------------------------------------------------
    # Step 1: Impute the data, Y:
    #----------------------------------------------------------------------------
    # Not implemented in this case
    #if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Not necessary here
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = tcrossprod(Y, t(Fmat)); sigma_tilde = sigma_et
    
    # Loop over each factor k = 1,...,K:
    for(k in 1:K){
      
      # Update the evoluation matrix:
      G_alpha[1,1] = ar_int[k]
      
      # Update the variances here:
      
      # Intercept/gamma:
      Wt[1,1,-T] = sigma_eta_tk[,k]^2; W0[1,1] = as.numeric(evolParams0_int$sigma_w0[k]^2)
      
      # Regression:
      if(p > 1){
        if(use_dynamic_reg){
          for(j in 1:(p-1)) Wt[-1, -1,][j,j,-T] = array(sigma_omega_tpk^2, c(T-1, p-1, K))[,j,k]
          diag(W0[-1, -1]) = as.numeric(matrix(evolParams0$sigma_w0^2, nrow = p-1)[,k])
        } else  W0[-1, -1] = diag(as.numeric(sigma_omega_kp[k,]^2))
      }
      
      # Sanity check for Wt: if variances too large, KFAS will stop running
      Wt[which(Wt > 10^6, arr.ind = TRUE)] = 10^6; W0[which(W0 > 10^6, arr.ind = TRUE)] = 10^6
      
      # Update the SSModel object given the new parameters
      kfas_model = update_kfas_model(Y.dlm = as.matrix(Y_tilde[,k] - mu_k[k]),
                                     Zt = X.arr,
                                     sigma_et = sigma_tilde,
                                     Gt = G_alpha, 
                                     Wt = Wt,
                                     W0 = W0,
                                     kfas_model = kfas_model)
      # Run the sampler
      alpha.arr[,,k] = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
      
      # Conditional mean from regression equation:
      Beta[,k] = mu_k[k] + rowSums(X*alpha.arr[,,k])
    }
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    # Not necessary, but update the conditional mean here
    Btheta = tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 5: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Btheta, svParams)
      sigma_et = svParams$sigma_et
    } else {
      # Or use uniform prior?
      sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Btheta)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, T)
    }
    #----------------------------------------------------------------------------
    # Step 6: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    # Centerend and non-centered:
    gamma_tk =  matrix(alpha.arr[,1,], nrow = T)
    
    gamma_tk_c = gamma_tk + mu_tk
    
    # Sample the unconditional mean term:
    mu_k = sampleARmu(yt = gamma_tk_c,
                      phi_j = ar_int,
                      sigma_tj = sigma_eta_tk,
                      priorPrec = 1/sigma_mu_k^2)
    mu_tk = matrix(rep(mu_k, each =  T), nrow = T)
    
    # And update the non-centered parameter:
    gamma_tk = gamma_tk_c - mu_tk
    
    # AR(1) coefficients:
    ar_int = sampleARphi(yt = gamma_tk,
                         phi_j = ar_int,
                         sigma_tj = sigma_eta_tk,
                         prior_phi = c(5,2))
    #prior_phi = NULL)
    
    # Then subtract the AR(1) part:
    eta_tk = gamma_tk[-1,] -  t(ar_int*t(gamma_tk[-T,]))
    
    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    a1_mu = uni.slice(a1_mu, g = function(a){
      dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_mu = uni.slice(a2_mu,g = function(a){
      sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                            delta.h = delta_eta_k,
                            a1 = a1_eta, a2 = a2_eta)
    sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
    a1_eta = uni.slice(a1_eta, g = function(a){
      dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_eta = uni.slice(a2_eta, g = function(a){
      sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Sample the corresponding prior variance term(s):
    evolParams_int = sampleEvolParams(omega = eta_tk/rep(sigma_delta_k, each = T-1), evolParams =  evolParams_int, evol_error = dist_reg_error)
    evolParams0_int = sampleEvol0(mu0 = as.matrix(gamma_tk[1,]), evolParams0 = evolParams0_int)
    
    # Update the error SD for gamma:
    sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T-1)
    
    # Cap at machine epsilon:
    sigma_eta_tk[which(sigma_eta_tk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
    #----------------------------------------------------------------------------
    # Step 7: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      
      if(use_dynamic_reg){
        
        # Dynamic setting
        
        # Regression (non-intercept) coefficients
        alpha_reg = matrix(alpha.arr[,-1,], nrow = T)
        
        # Random walk, so compute difference for innovations:
        omega = diff(alpha_reg)
        
        #----------------------------------------------------------------------------
        if(add_p_shrink){
          # First part: j-specific shrinkage
          # Scale omega by lambda_omega_tpk:
          omega_scale = omega/lambda_omega_tpk
          # Sum of squares (over t,k)
          ss_omega_j = apply(array(omega_scale^2, c(T-1, p-1, K)), 2, sum)
          # Offset:
          ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
          
          # SD term:
          lambda_j = 1/sqrt(rgamma(n = p - 1,
                                   shape = 1/2 + (T-1)*K/2,
                                   rate = px_lambda_j + ss_omega_j/2))
          
          # Sample the PX, global parameters outside this if-statement
          
          # Recurring term:
          rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K)
        }
        #----------------------------------------------------------------------------
        # Second part: tjk-specific shrinkage
        
        # Scale by lambda_j
        omega_scale = omega/rep_lambda_j
        
        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)
        
        # SD term:
        lambda_omega_tpk = evolParams$sigma_wt
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters
        
        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j
        
        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
      } else{
        
        # Non-dynamic setting: grab the first one (all the same),
        # and store as (K x p-1) matrix
        #omega = t(alpha.arr[1, -1, ])
        omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)
        
        #----------------------------------------------------------------------------
        if(add_p_shrink){
          # First part: j-specific shrinkage
          
          # Sum of squares (over k)
          ss_omega_j = colSums((omega/lambda_omega_kp)^2)
          # Offset:
          ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
          
          # SD term:
          lambda_j = 1/sqrt(rgamma(n = p - 1,
                                   shape = 1/2 + K/2,
                                   rate = px_lambda_j + ss_omega_j/2))
          # Recurring term:
          rep_lambda_j = rep(lambda_j, each = K)
        }
        #----------------------------------------------------------------------------
        # Second part: jk-specific shrinkage
        
        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega/rep_lambda_j, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        
        # SD term:
        lambda_omega_kp = evolParams$sigma_wt
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters
        
        # SD for omega:
        sigma_omega_kp = lambda_omega_kp*rep_lambda_j
        
        # Cap at machine epsilon:
        sigma_omega_kp[which(sigma_omega_kp < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
      }
      #----------------------------------------------------------------------------
      if(add_p_shrink){
        # In both dynamic/non-dynamic cases, need to sample the PX-parameters
        # for lambda_j and global shrinkage parameters
        
        # PX-j:
        px_lambda_j = rgamma(n = p - 1, shape = 1, rate = 1/lambda_j^2 + 1/lambda_0^2)
        
        # Global shrinkage:
        # SD term:
        lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_lambda_0 + sum(px_lambda_j)))
        
        # PX:
        px_lambda_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    # Not necessary here
    #if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}
    
    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,,] = alpha.arr
        if(!is.na(match('mu_k', mcmc_params))) post.mu_k[isave,] = mu_k
        if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et[isave,] = sigma_et
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Btheta # + sigma_e*rnorm(length(Y))
        if(!is.na(match('ar_phi', mcmc_params))) post.ar_phi[isave,] = ar_int
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Btheta), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }
  
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('mu_k', mcmc_params))) mcmc_output$mu_k = post.mu_k
  if(!is.na(match('sigma_et', mcmc_params))) mcmc_output$sigma_et = post.sigma_et
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('ar_phi', mcmc_params))) mcmc_output$ar_phi = post.ar_phi
  
  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_et), m),
                            log = TRUE), na.rm=TRUE)
    
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d
    
    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))
  
  return (mcmc_output);
}
#' MCMC Sampling Algorithm for the Dynamic Function-on-Scalars Regression Model with Random Walk Errors
#'
#' Runs the MCMC for the dynamic function-on-scalars regression model based on
#' an FDLM-type expansion. Here, we assume the factor regression has RW errors.
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param use_fpca logical; if TRUE, use FPCA, otherwise use orthonormal B-splines
#' @param use_dynamic_reg logical; if TRUE, regression coefficients are dynamic
#' (with random walk models), otherwise independent
#' @param dist_reg_coef prior distribution for the regression coefficient evolution errors;
#' must be one of
#' \itemize{
#' \item "NIG" (Normal-inverse-Gamma)
#' \item "HS" (horseshoe prior)
#' \item "DHS" (dynamic horseshoe prior)
#' \item "BL" (Bayesian lasso)
#' \item "SV" (Gaussian stochastic volatility)
#' }
#' @param dist_reg_error prior distribution for the regression-level errors; same options as \code{dist_reg_coef}
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "alpha" (regression coefficients)
#' \item "mu_k" (intercept term for factor k)
#' \item "sigma_et" (observation error SD; possibly dynamic)
#' \item "Yhat" (fitted values)
#' }
#' @param add_p_shrink logical; when TRUE, include a predictor-specific shrinkage term
#' (i.e., the scale term will be half-Cauchy)
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model
#' for the observation error variance
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{Tm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x T x m},  may be inefficient
#'
#' @examples
#' # Fixme
#'
#' @import dsp KFAS truncdist
#' @export
dfosr_rw_basis = function(Y, tau, X = NULL, use_fpca = FALSE,
                    use_dynamic_reg = TRUE,
                    dist_reg_coef = "DHS", dist_reg_error = "NIG",
                    nsave = 1000, nburn = 1000, nskip = 10,
                    mcmc_params = list("beta", "fk", "alpha"),
                    add_p_shrink = FALSE,
                    use_obs_SV = FALSE,
                    computeDIC = TRUE){
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)
  
  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y); d = ncol(tau)
  
  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  #----------------------------------------------------------------------------
  # Initialize the main terms
  
  if(any(is.na(Y))) stop("Missing data not implemented for basis methods")
  Yna = Y # for consistency later  
  
  if(use_fpca){
    if(d > 1) stop("FPCA only implemented for d = 1")
    
    Fmat = fpca.face(Y, center = TRUE, 
                     argvals = as.numeric(tau01),
                     knots=  max(20, min(ceiling(floor(median(rowSums(!is.na(Y))))/4), 150)),
                     pve=0.99)$efunctions
  } else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  #} else Fmat = getSplineInfo_d(tau = tau01, m_eff = floor(0.1*median(rowSums(!is.na(Y)))), orthonormalize = TRUE)$Bmat
  
  # Initialize the factors
  Beta = tcrossprod(Y, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value
  
  # Initialize the conditional expectation:
  Btheta = tcrossprod(Beta, Fmat)
  
  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }
  #----------------------------------------------------------------------------
  # Predictors:
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, T), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors:
  p = ncol(X)
  
  # This array is useful for setting up the DLM
  X.arr = array(0, c(K, K*p, T)); for(i in 1:T) X.arr[,,i] = blockDiag(matrix(X[i,], nrow=1) ,K)
  
  # Initialize the SSModel:
  kfas_model = update_kfas_model(Y.dlm = Beta, Zt = X.arr)
  
  # Indices identifying the intercept (or gamma_k term):
  ind.intercept = seq(1, K*p, by = p)
  
  # Identify the dynamic/non-dynamic components:
  if(!use_dynamic_reg) diag(kfas_model$R[,,1])[-ind.intercept] = 0
  #----------------------------------------------------------------------------
  # Initialize the regression terms:
  
  # Evolution matrix (must be identity for static components)
  G_alpha = diag(K*p)
  
  # Update the SSModel object given the new parameters
  kfas_model = update_kfas_model(Y.dlm = Beta, Zt = X.arr,Gt = G_alpha, kfas_model = kfas_model)
  
  # Run the sampler
  alpha = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
  
  # Store as an array:
  alpha.arr = array(alpha, c(T, p, K))
  
  # Conditional mean from regression equation:
  for(k in 1:K) Beta[,k] = rowSums(X*alpha.arr[,,k])
  #----------------------------------------------------------------------------
  # Evolution error variance:
  Wt = array(diag(K*p), c(K*p, K*p, T)); W0 = diag(10^-4, K*p);
  
  # Intercept (or gamma) components:
  gamma_tk =  as.matrix(alpha[,ind.intercept])
  
  # Difference:
  eta_tk = diff(gamma_tk)
  
  # Initialize the corresponding prior variance term(s):
  evolParams_int = initEvolParams(omega = eta_tk, evol_error = dist_reg_error)
  evolParams0_int = initEvol0(mu0 = as.matrix(gamma_tk[1,]))
  
  # MGP term:
  a1_eta = 2; a2_eta = 3
  delta_eta_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
  
  # Update the error SD for gamma:
  sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T-1)
  
  # Evolution and initial variance:
  for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_eta_tk[,k]^2
  diag(W0[ind.intercept, ind.intercept]) = as.numeric(evolParams0_int$sigma_w0^2)
  #----------------------------------------------------------------------------
  # Non-intercept term:
  if(p > 1){
    
    if(use_dynamic_reg){
      
      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])
      
      # Innovation:
      omega = diff(alpha_reg)
      
      # First, initialize the variable selection scaling terms:
      if(add_p_shrink){
        lambda_j = apply(array(omega, c(T-1, p-1, K)), 2, sd)
        lambda_0 = median(lambda_j) # Global term
        px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
        rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K) # Recurring term
      } else rep_lambda_j = rep(1, (T-1)*(p-1)*K) # Recurring scale term: set to 1
      
      # Initialize the variance components:
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))
      # This is the scaling SD
      lambda_omega_tpk = evolParams$sigma_wt
      
      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j
      
      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)
      
    } else{
      
      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      omega = t(alpha.arr[1, -1, ])
      
      # First, initialize the variable selection scaling terms:
      if(add_p_shrink){
        lambda_j = apply(omega, 2, sd)
        lambda_0 = median(lambda_j) # Global term
        px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
        rep_lambda_j = rep(lambda_j, each = K) # Recurring term
      } else rep_lambda_j = rep(1, (p-1)*K) # Recurring scale term
      
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      
      # This is the scaling SD (Note: doesn't actually depend on T)
      lambda_omega_pk = evolParams$sigma_wt
      
      # SD for omega:
      sigma_omega_pk = lambda_omega_pk*rep_lambda_j
      
      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))
    }
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, T, p, K))
  if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et = array(NA, c(nsave, T))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, T, m))
  if(computeDIC) post_loglike = numeric(nsave)
  
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    #----------------------------------------------------------------------------
    # Step 1: Impute the data, Y:
    #----------------------------------------------------------------------------
    # Not implemented in this case
    #if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Not necessary here
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = tcrossprod(Y, t(Fmat)); sigma_tilde = sigma_et
    
    # Sanity check for Wt: if variances too large, KFAS will stop running
    Wt[which(Wt > 10^6, arr.ind = TRUE)] = 10^6; W0[which(W0 > 10^6, arr.ind = TRUE)] = 10^6
    
    # Update the SSModel object and sample:
    kfas_model = update_kfas_model(Y.dlm = Y_tilde,
                                   Zt = X.arr,
                                   sigma_et = sigma_tilde,
                                   Gt = G_alpha,
                                   Wt = Wt,
                                   W0 = W0,
                                   kfas_model = kfas_model)
    alpha = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
    
    # Store as an array:
    alpha.arr = array(alpha, c(T, p, K))
    
    # Update the factors:
    for(k in 1:K) Beta[,k] = rowSums(X*alpha.arr[,,k])
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    # Not necessary, but update the conditional mean here
    Btheta = tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 5: Sample the observation error variance
    #----------------------------------------------------------------------------
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Btheta, svParams)
      sigma_et = svParams$sigma_et
    } else {
      # Or use uniform prior?
      sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Btheta)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, T)
    }
    #----------------------------------------------------------------------------
    # Step 6: Sample the intercept/gamma parameters
    #----------------------------------------------------------------------------
    gamma_tk =  as.matrix(alpha[,ind.intercept])
    
    # Difference:
    eta_tk = diff(gamma_tk)
    
    # Prior variance: MGP
    # Standardize, then reconstruct as matrix of size T x K:
    delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                            delta.h = delta_eta_k,
                            a1 = a1_eta, a2 = a2_eta)
    sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
    a1_eta = uni.slice(a1_eta, g = function(a){
      dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
    a2_eta = uni.slice(a2_eta, g = function(a){
      sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    
    # Sample the corresponding prior variance term(s):
    evolParams_int = sampleEvolParams(omega = eta_tk/rep(sigma_delta_k, each = T-1), evolParams =  evolParams_int, evol_error = dist_reg_error)
    evolParams0_int = sampleEvol0(mu0 = as.matrix(gamma_tk[1,]), evolParams0 = evolParams0_int)
    
    # Update the error SD for gamma:
    sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T-1)
    
    # Cap at machine epsilon:
    sigma_eta_tk[which(sigma_eta_tk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
    
    # Evolution and initial variance:
    for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_eta_tk[,k]^2
    diag(W0[ind.intercept, ind.intercept]) = as.numeric(evolParams0_int$sigma_w0^2)
    #----------------------------------------------------------------------------
    # Step 7: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      
      if(use_dynamic_reg){
        
        # Dynamic setting
        
        # Regression (non-intercept) coefficients
        alpha_reg = as.matrix(alpha[,-ind.intercept])
        
        # Random walk, so compute difference for innovations:
        omega = diff(alpha_reg)
        
        #----------------------------------------------------------------------------
        if(add_p_shrink){
          # First part: j-specific shrinkage
          # Scale omega by lambda_omega_tpk:
          omega_scale = omega/lambda_omega_tpk
          # Sum of squares (over t,k)
          ss_omega_j = apply(array(omega_scale^2, c(T-1, p-1, K)), 2, sum)
          # Offset:
          ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
          
          # SD term:
          lambda_j = 1/sqrt(rgamma(n = p - 1,
                                   shape = 1/2 + (T-1)*K/2,
                                   rate = px_lambda_j + ss_omega_j/2))
          
          # Sample the PX, global parameters outside this if-statement
          
          # Recurring term:
          rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K)
        }
        #----------------------------------------------------------------------------
        # Second part: tjk-specific shrinkage
        
        # Scale by lambda_j
        omega_scale = omega/rep_lambda_j
        
        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)
        
        # SD term:
        lambda_omega_tpk = evolParams$sigma_wt
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters
        
        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j
        
        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
        
        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)
        
      } else{
        
        # Non-dynamic setting: grab the first one (all the same),
        # and store as (K x p-1) matrix
        omega = t(alpha.arr[1, -1, ])
        
        #----------------------------------------------------------------------------
        if(add_p_shrink){
          # First part: j-specific shrinkage
          
          # Sum of squares (over k)
          ss_omega_j = colSums((omega/lambda_omega_pk)^2)
          # Offset:
          ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8
          
          # SD term:
          lambda_j = 1/sqrt(rgamma(n = p - 1,
                                   shape = 1/2 + K/2,
                                   rate = px_lambda_j + ss_omega_j/2))
          # Recurring term:
          rep_lambda_j = rep(lambda_j, each = K)
        }
        #----------------------------------------------------------------------------
        # Second part: jk-specific shrinkage
        
        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega/rep_lambda_j, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        
        # SD term:
        lambda_omega_pk = evolParams$sigma_wt
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters
        
        # SD for omega:
        sigma_omega_pk = lambda_omega_pk*rep_lambda_j
        
        # Cap at machine epsilon:
        sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
        
        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))
      }
      #----------------------------------------------------------------------------
      if(add_p_shrink){
        # In both dynamic/non-dynamic cases, need to sample the PX-parameters
        # for lambda_j and global shrinkage parameters
        
        # PX-j:
        px_lambda_j = rgamma(n = p - 1, shape = 1, rate = 1/lambda_j^2 + 1/lambda_0^2)
        
        # Global shrinkage:
        # SD term:
        lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_lambda_0 + sum(px_lambda_j)))
        
        # PX:
        px_lambda_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    # Not necessary
    #if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}
    
    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,,] = alpha.arr
        if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et[isave,] = sigma_et
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Btheta # + sigma_e*rnorm(length(Y))
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Btheta), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }
  
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('sigma_et', mcmc_params))) mcmc_output$sigma_et = post.sigma_et
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  
  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_et), m),
                            log = TRUE), na.rm=TRUE)
    
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d
    
    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))
  
  return (mcmc_output);
}
