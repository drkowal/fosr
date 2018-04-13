# Easy way to set the parameters:
#if(FALSE){
#  library(KFAS)
#  K = 6; use_dynamic_reg = TRUE; dist_reg_coef = "DHS"; dist_reg_error = "NIG";
#  nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi", "Yhat")
#  use_obs_SV = FALSE; includeBasisInnovation = FALSE; computeDIC = TRUE
#  add_p_shrink = TRUE
#
#  #use_dynamic_reg = FALSE
#  #for(j in 1:p) plot(as.zoo(alpha.arr[,j,]), ylim = range(alpha.arr), main = j)
#}

naive_lm = function(Y, X = NULL,
                    nsave = 1000, nburn = 1000, nskip = 3,
                    mcmc_params = list("alpha", "Yhat"),
                    computeDIC = TRUE){

  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y);

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

  # p x m regression "functions"
  alpha_pm = matrix(0, nrow = p, ncol = m)
  for(j in 1:m) alpha_pm[,j] = lm(Y[,j] ~ X - 1)$coef

  # SD of regression "functions"
  sigma_alpha = sd(alpha_pm[-1,])

  # SD of intercept:
  sigma_mu = sd(alpha_pm[1,])

  # Fitted values:
  Yhat = tcrossprod(X, t(alpha_pm))

  # SD of error:
  sigma_e = sd(Y - Yhat)

  # Some useful terms:
  XtX = crossprod(X);
  XtY = crossprod(X,Y)

  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, m))
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
    # Step 1: Sample the regression coefficients
    #----------------------------------------------------------------------------
    chQ = chol(XtX/sigma_e^2 + diag(1/c(sigma_mu, rep(sigma_alpha, p-1))^2))
    ell_j = XtY/sigma_e^2
    #for(j in 1:m) alpha_pm[,j] = backsolve(chQ, forwardsolve(t(chQ), ell_j[,j]) + rnorm(p))
    alpha_pm = backsolve(chQ,forwardsolve(t(chQ), ell_j) + rnorm(length(ell_j)))

    # Fitted values:
    Yhat = tcrossprod(X, t(alpha_pm))
    #----------------------------------------------------------------------------
    # Step 2: Sample the variances
    #----------------------------------------------------------------------------
    sigma_e = 1/sqrt(rgamma(n = 1, shape = T*m/2, rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_alpha = 1/sqrt(rgamma(n = 1, shape = 0.01 + (p-1)*m/2, rate = 0.01 + sum(alpha_pm[-1,]^2)/2))
    sigma_mu = 1/sqrt(rgamma(n = 1, shape = 0.01 + m/2, rate = 0.01 + sum(alpha_pm[1,]^2)/2))

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = alpha_pm
        if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et[isave,] = rep(sigma_e, T)
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Yhat
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Y), mean = matrix(Yhat), sd = rep(sigma_e, T*m), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Y),
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

  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
#' MCMC Sampling Algorithm for the Function-on-Scalars Regression Model
#'
#' Runs the MCMC for the function-on-scalars regression model based on
#' an FDLM-type expansion. Here we assume the factor regression has independent errors.
#' In addition, the prior on the regression coefficients is a multiplicative gamma prior (MGP).
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{T x p} matrix of predictors; if NULL, only include an intercept
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
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
#' \item "sigma_et" (observation error SD)
#' \item "Yhat" (fitted values)
#' }
#' @param includeBasisInnovation logical; when TRUE, include an iid basis coefficient term for residual correlation
#' (i.e., the idiosyncratic error term for a factor model on the full basis matrix)
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
fosr_mgp = function(Y, tau, X = NULL, K = NULL,
                dist_reg_error = "NIG",
                nsave = 1000, nburn = 1000, nskip = 3,
                mcmc_params = list("beta", "fk", "alpha"),
                includeBasisInnovation = FALSE,
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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  K = ncol(Beta) # to be sure we have the right value

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    nu = t(BtY) - BetaPsit; sigma_nu = sd(nu)
    theta = BetaPsit + nu; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_nu = 0

  # Initialize the (time-dependent) observation error SD:
  sigma_e = sd(Y - Btheta, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
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

  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:K) {
    if(p >= T){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_et,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Y, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_tk[,k] = Beta[,k] - X%*%alpha_pk[,k]
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
    # MGP term for each predictor:
    a1_omega = 2; a2_omega = 3
    delta_omega_pk = tau_omega_pk = sigma_omega_pk = matrix(1, nrow = p-1, ncol = K)

    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

    sigma_omega_pk = abs(omega)
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2 + sigma_nu^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    if(includeBasisInnovation){
      Y_tilde =  tcrossprod(theta, t(Psi)); sigma_tilde = rep(sigma_nu, T)
    } else {
      Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et
    }

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 + sigma_gamma_tk[,k]^2)

      if(p >= T){
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        ell_k = crossprod(X, y_tilde_k/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }

    # And sample the errors gamma_tk:
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + matrix(1/sigma_gamma_tk^2))
    postMean = matrix((Y_tilde - X%*%alpha_pk)/rep(sigma_tilde^2, times = K))*postSD^2
    gamma_tk = matrix(rnorm(n = T*K, mean = postMean, sd = postSD), nrow = T)

    # Update the factors:
    Beta = X%*%alpha_pk + gamma_tk

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      chQtheta = sqrt(sigma_e^-2 + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
      linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_nu^2 # Linear term from the posterior
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)

      sigma_nu = 1/sqrt(truncdist::rtrunc(1, "gamma",
                                          a = 10^-8, b = Inf,
                                          shape = (length(theta)+1)/2,
                                          rate = 1/2*sum((theta - BetaPsit)^2)))
    } else {theta = BetaPsit; sigma_nu = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}
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
      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      #----------------------------------------------------------------------------
      # jk-specific shrinkage: loop over each predictor
      for(j in 1:(p-1)){
        delta_omega_pk[j,] = sampleMGP(theta.jh = omega[j,],
                                       delta.h = delta_omega_pk[j,],
                                       a1 = a1_omega, a2 = a2_omega)
        tau_omega_pk[j,] = cumprod(delta_omega_pk[j,])
      }
      # Update the MGP parameters, common to all j:
      a1_omega = uni.slice(a1_omega, g = function(a){
        sum(dgamma(delta_omega_pk[,1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_omega = uni.slice(a2_omega, g = function(a){
        sum(dgamma(delta_omega_pk[,-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
      #----------------------------------------------------------------------------
      # Last part: update the variance parameters

      # SD for omega:
      sigma_omega_pk = 1/(sqrt(tau_omega_pk))

      # Cap at machine epsilon:
      sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)
     }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

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
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
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

  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
# dfosr_ar w/ NEW options:
  # no MGP: j-specific shrinkage + (jkt)-specific shrinkage
  # MGP: use MGP prior for each predictor, j, + (jkt)-specific shrinkage
# The point is that MGP may not be sensible for predictors:
  # what if only one k = k_* is useful for predictor x_j?

# Modified dfosr_ar for HS priors:
  # lambda_(jkt)^-2 ~ Ga(3/2, 3/2)
  # tau_(jk) ~ MGP
dfosr_ar_mgp = function(Y, tau, X = NULL, K = NULL,
                       use_dynamic_reg = TRUE,
                       dist_reg_error = "NIG",
                       nsave = 1000, nburn = 1000, nskip = 10,
                       mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                       use_obs_SV = FALSE,
                       includeBasisInnovation = FALSE,
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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    nu = t(BtY) - BetaPsit; sigma_nu = sd(nu)
    theta = BetaPsit + nu; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_nu = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
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

    # MGP term for each predictor:
    a1_omega = 2; a2_omega = 3
    delta_omega_pk = tau_omega_pk = matrix(1, nrow = p-1, ncol = K)

    if(use_dynamic_reg){

      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])

      # Innovation:
      omega = diff(alpha_reg)

      # Initialize the (jkt)-specific component scaling SD:
      lambda_omega_tpk = matrix(1, nrow = (T-1), ncol = (p-1)*K)

      # Initialize the initial variance components:
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))

      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep(matrix(1/sqrt(tau_omega_pk)), each = (T-1))

      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

    } else{

      # Non-dynamic setting: grab the first one (all the same)
      omega = alpha.arr[1, -1, ]

      sigma_omega_pk = 1/(sqrt(tau_omega_pk))

      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(sigma_omega_pk^2))
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2 + sigma_nu^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    if(includeBasisInnovation){
      Y_tilde =  tcrossprod(theta, t(Psi)); sigma_tilde = sigma_nu
    } else {
      Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et
    }

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

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_nu^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_nu^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)

      sigma_nu = 1/sqrt(truncdist::rtrunc(1, "gamma",
                                          a = 10^-8, b = Inf,
                                          shape = (length(theta)+1)/2,
                                          rate = 1/2*sum((theta - BetaPsit)^2)))
    } else {theta = BetaPsit; sigma_nu = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}
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
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
    a2_mu = uni.slice(a2_mu,g = function(a){
      sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                            delta.h = delta_eta_k,
                            a1 = a1_eta, a2 = a2_eta)
    sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
    a1_eta = uni.slice(a1_eta, g = function(a){
      dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
    a2_eta = uni.slice(a2_eta, g = function(a){
      sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

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

        # Initial variance sampler:
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)

        #----------------------------------------------------------------------------
        # First part: tjk-specific shrinkage

        # Scale by tau_omega_pk:
        omega_scale = omega*rep(matrix(sqrt(tau_omega_pk)), each = (T-1))

        # Squared, then check for numerical issues:
        ss_omega_scale = omega_scale^2; ss_omega_scale = ss_omega_scale + (ss_omega_scale < 10^-16)*10^-8

        # Sample the SD piece:
        nu = 3
        lambda_omega_tpk = 1/sqrt(matrix(rgamma(n = (T-1)*(p-1)*K,
                                                shape = (nu + 1)/2,
                                                rate = nu/2 + ss_omega_scale/2), nrow = T-1))
        #----------------------------------------------------------------------------
        # Second part: (jk)-specific shrinkage (MGP)
        omega_scale = array(omega/lambda_omega_tpk, c(T-1, p-1, K))

        # For each predictor:
        for(j in 1:(p-1)){
          delta_omega_pk[j,] = sampleMGP(theta.jh = omega_scale[,j,],
                                         delta.h = delta_omega_pk[j,],
                                         a1 = a1_omega, a2 = a2_omega)
          tau_omega_pk[j,] = cumprod(delta_omega_pk[j,])
        }
        # Update the MGP parameters, common to all j:
        a1_omega = uni.slice(a1_omega, g = function(a){
          sum(dgamma(delta_omega_pk[,1], shape = a, rate = 1, log = TRUE)) +
            dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
        a2_omega = uni.slice(a2_omega, g = function(a){
          sum(dgamma(delta_omega_pk[,-1], shape = a, rate = 1, log = TRUE)) +
            dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters

        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep(matrix(1/sqrt(tau_omega_pk)), each = (T-1))

        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      } else{

        # Non-dynamic setting: grab the first one (all the same)
        omega = alpha.arr[1, -1, ]
        #----------------------------------------------------------------------------
        # jk-specific shrinkage: loop over each predictor
        for(j in 1:(p-1)){
          delta_omega_pk[j,] = sampleMGP(theta.jh = omega[j,],
                                         delta.h = delta_omega_pk[j,],
                                         a1 = a1_omega, a2 = a2_omega)
          tau_omega_pk[j,] = cumprod(delta_omega_pk[j,])
        }
        # Update the MGP parameters, common to all j:
        a1_omega = uni.slice(a1_omega, g = function(a){
          sum(dgamma(delta_omega_pk[,1], shape = a, rate = 1, log = TRUE)) +
            dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
        a2_omega = uni.slice(a2_omega, g = function(a){
          sum(dgamma(delta_omega_pk[,-1], shape = a, rate = 1, log = TRUE)) +
            dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
        #----------------------------------------------------------------------------
        # Last part: update the variance parameters

        # SD for omega:
        sigma_omega_pk = 1/(sqrt(tau_omega_pk))

        # Cap at machine epsilon:
        sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(sigma_omega_pk^2))
        #----------------------------------------------------------------------------
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

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

# Modified dfosr_ar for HS priors:
# lambda_(jkt)^-2 ~ Ga(3/2, 3/2)
# tau_(jk) ~ MGP+HS
dfosr_ar_mgp_hs = function(Y, tau, X = NULL, K = NULL,
                       use_dynamic_reg = TRUE,
                       dist_reg_error = "NIG",
                       nsave = 1000, nburn = 1000, nskip = 10,
                       mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                       use_obs_SV = FALSE,
                       includeBasisInnovation = FALSE,
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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    nu = t(BtY) - BetaPsit; sigma_nu = sd(nu)
    theta = BetaPsit + nu; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_nu = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
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

    # MGP term for each predictor:
    a_omega = 3 # Only need this piece
    delta_omega_pk = tau_omega_pk = matrix(1, nrow = p-1, ncol = K)
    px_omega_p = rep(1, p-1) # Parameter-expansion variable for each p

    # Global scale parameter + px parameter:
    lambda_0 = px_omega_0 = 1

    if(use_dynamic_reg){

      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])

      # Innovation:
      omega = diff(alpha_reg)

      # Initialize the (jkt)-specific component scaling SD:
      lambda_omega_tpk = matrix(1, nrow = (T-1), ncol = (p-1)*K)

      # Initialize the initial variance components:
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))

      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep(matrix(1/sqrt(tau_omega_pk)), each = (T-1))

      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

    } else{

      # Non-dynamic setting: grab the first one (all the same)
      omega = alpha.arr[1, -1, ]

      sigma_omega_pk = 1/(sqrt(tau_omega_pk))

      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(sigma_omega_pk^2))
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2 + sigma_nu^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    if(includeBasisInnovation){
      Y_tilde =  tcrossprod(theta, t(Psi)); sigma_tilde = sigma_nu
    } else {
      Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et
    }

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

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_nu^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_nu^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)

      sigma_nu = 1/sqrt(truncdist::rtrunc(1, "gamma",
                                          a = 10^-8, b = Inf,
                                          shape = (length(theta)+1)/2,
                                          rate = 1/2*sum((theta - BetaPsit)^2)))
    } else {theta = BetaPsit; sigma_nu = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}
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
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
    a2_mu = uni.slice(a2_mu,g = function(a){
      sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                            delta.h = delta_eta_k,
                            a1 = a1_eta, a2 = a2_eta)
    sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
    a1_eta = uni.slice(a1_eta, g = function(a){
      dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
    a2_eta = uni.slice(a2_eta, g = function(a){
      sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

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

        # Initial variance sampler:
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)

        #----------------------------------------------------------------------------
        # First part: tjk-specific shrinkage

        # Scale by tau_omega_pk:
        omega_scale = omega*rep(matrix(sqrt(tau_omega_pk)), each = (T-1))

        # Squared, then check for numerical issues:
        ss_omega_scale = omega_scale^2; ss_omega_scale = ss_omega_scale + (ss_omega_scale < 10^-16)*10^-8

        # Sample the SD piece:
        nu = 3
        lambda_omega_tpk = 1/sqrt(matrix(rgamma(n = (T-1)*(p-1)*K,
                                                shape = (nu + 1)/2,
                                                rate = nu/2 + ss_omega_scale/2), nrow = T-1))
        #----------------------------------------------------------------------------
        # Second part: (jk)-specific shrinkage (MGP)
        omega_scale = array(omega/lambda_omega_tpk, c(T-1, p-1, K))

        # For each predictor:
        for(j in 1:(p-1)){
          # SS scaled omega (for each k, summed across t for predictor j)
          ss_omega_k = colSums(omega_scale[,j,]^2)

          # k = 1 case is separate:
          tau_not_1 = cumprod(delta_omega_pk[j,])/delta_omega_pk[j,1]
          ss_omega_1 = sum(tau_not_1*ss_omega_k); ss_omega_1 = ss_omega_1 + (ss_omega_1 < 10^-16)*10^-8
          delta_omega_pk[j,1] = rgamma(n = 1,
                                       shape = 1/2 + (T-1)*K/2,
                                       rate = px_omega_p[j] + ss_omega_1/2)
          # k > 1:
          if(K > 1){for(h in 2:K){
            tau_not_h = cumprod(delta_omega_pk[j,])/delta_omega_pk[j,h]
            ss_omega_h = sum(tau_not_h[h:K]*ss_omega_k[h:K]); ss_omega_h = ss_omega_h + (ss_omega_h < 10^-16)*10^-8
            delta_omega_pk[j,h] = rgamma(n = 1,
                                         shape = a_omega + (T-1)/2*(K - h + 1),
                                         rate = 1 + 1/2*ss_omega_h)
          }}
          tau_omega_pk[j,] = cumprod(delta_omega_pk[j,])
        }
        # Sample the PX components:
        px_omega_p = rgamma(n = p-1, shape = 1, rate = delta_omega_pk[,1] + 1/lambda_0^2)

        # And the global parameters (w/ PX):
        lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_omega_0 + sum(px_omega_p)))
        px_omega_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)

        # Shape parameter for k > 2:
        a_omega = uni.slice(a_omega, g = function(a){
          sum(dgamma(delta_omega_pk[,-1], shape = a, rate = 1, log = TRUE)) +
            dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters

        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep(matrix(1/sqrt(tau_omega_pk)), each = (T-1))

        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      } else{

        # Non-dynamic setting: grab the first one (all the same)
        omega = alpha.arr[1, -1, ]
        #----------------------------------------------------------------------------
        # jk-specific shrinkage: loop over each predictor
        for(j in 1:(p-1)){
          # SS scaled omega (for each k, summed across t for predictor j)
          ss_omega_k = omega[j,]^2

          # k = 1 case is separate:
          tau_not_1 = cumprod(delta_omega_pk[j,])/delta_omega_pk[j,1]
          ss_omega_1 = sum(tau_not_1*ss_omega_k); ss_omega_1 = ss_omega_1 + (ss_omega_1 < 10^-16)*10^-8
          delta_omega_pk[j,1] = rgamma(n = 1,
                                       shape = 1/2 + K/2,
                                       rate = px_omega_p[j] + ss_omega_1/2)
          # k > 1:
          if(K > 1){for(h in 2:K){
            tau_not_h = cumprod(delta_omega_pk[j,])/delta_omega_pk[j,h]
            ss_omega_h = sum(tau_not_h[h:K]*ss_omega_k[h:K]); ss_omega_h = ss_omega_h + (ss_omega_h < 10^-16)*10^-8
            delta_omega_pk[j,h] = rgamma(n = 1,
                                         shape = a_omega + (K - h + 1)/2,
                                         rate = 1 + 1/2*ss_omega_h)
          }}
          tau_omega_pk[j,] = cumprod(delta_omega_pk[j,])
        }
        # Sample the PX components:
        px_omega_p = rgamma(n = p-1, shape = 1, rate = delta_omega_pk[,1] + 1/lambda_0^2)

        # And the global parameters (w/ PX):
        lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_omega_0 + sum(px_omega_p)))
        px_omega_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)

        # Shape parameter for k > 2:
        a_omega = uni.slice(a_omega, g = function(a){
          sum(dgamma(delta_omega_pk[,-1], shape = a, rate = 1, log = TRUE)) +
            dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
        #----------------------------------------------------------------------------
        # Last part: update the variance parameters

        # SD for omega:
        sigma_omega_pk = 1/(sqrt(tau_omega_pk))

        # Cap at machine epsilon:
        sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(sigma_omega_pk^2))
        #----------------------------------------------------------------------------
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

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
# Modified dfosr_ar for HS priors:
  # lambda_(jkt) ~ C+(0,1)
  # tau_(jk)^(-1/2) ~ C+(0, lambda_j)
  # lambda_j ~ C+(0, lambda_0)
  # lambda_0 ~ C+(0, 1)
dfosr_ar_hs = function(Y, tau, X = NULL, K = NULL,
                       use_dynamic_reg = TRUE,
                       dist_reg_error = "NIG",
                       nsave = 1000, nburn = 1000, nskip = 10,
                       mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                       use_obs_SV = FALSE,
                       includeBasisInnovation = FALSE,
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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    nu = t(BtY) - BetaPsit; sigma_nu = sd(nu)
    theta = BetaPsit + nu; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_nu = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
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

      # Hierarchical half-Cauchy priors

      # tpk component:
      lambda_omega_tpk = abs(omega);
      px_lambda_tpk = matrix(1, nrow = (T-1), ncol = (p-1)*K)

      # pk components:
      lambda_omega_pk = matrix(colMeans(lambda_omega_tpk), nrow = p-1, ncol = K);
      px_lambda_pk = matrix(1, nrow = (p-1), ncol =K)

      # p components:
      lambda_omega_p = rowMeans(lambda_omega_pk)
      px_lambda_p = rep(1, p-1)

      # Global components:
      lambda_omega_0 = mean(lambda_omega_p)
      px_lambda_0 = 1

      # Initialize the initial variance components:
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))

      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep(matrix(lambda_omega_pk), each = (T-1))

      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

    } else{

      stop("non-dynamic case still in progress, and not really necessary")

      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      #omega = t(alpha.arr[1, -1, ])
      omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

      evolParams = initEvolParams(omega = omega, evol_error = dist_reg_coef)

      # This is the scaling SD (Note: doesn't actually depend on T)
      lambda_omega_kp = evolParams$sigma_wt

      # SD for omega:
      sigma_omega_pk = t(lambda_omega_kp)/sqrt(tau_omega_pk)

      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(sigma_omega_pk^2))
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2 + sigma_nu^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    if(includeBasisInnovation){
      Y_tilde =  tcrossprod(theta, t(Psi)); sigma_tilde = sigma_nu
    } else {
      Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et
    }

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

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_nu^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_nu^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)

      sigma_nu = 1/sqrt(truncdist::rtrunc(1, "gamma",
                                          a = 10^-8, b = Inf,
                                          shape = (length(theta)+1)/2,
                                          rate = 1/2*sum((theta - BetaPsit)^2)))
    } else {theta = BetaPsit; sigma_nu = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}
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
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
    a2_mu = uni.slice(a2_mu,g = function(a){
      sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

    # Variance part:
    # Standardize, then reconstruct as matrix of size T x K:
    delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                            delta.h = delta_eta_k,
                            a1 = a1_eta, a2 = a2_eta)
    sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
    a1_eta = uni.slice(a1_eta, g = function(a){
      dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
    a2_eta = uni.slice(a2_eta, g = function(a){
      sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

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

        # Initial variance sampler:
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)

        #----------------------------------------------------------------------------
        # First part: tjk-specific shrinkage

        # Scale by lambda_omega_pk
        omega_scale = omega/rep(matrix(lambda_omega_pk), each = (T-1))

        # Squared, then check for numerical issues:
        ss_omega_scale = omega_scale^2; ss_omega_scale = ss_omega_scale + (ss_omega_scale < 10^-16)*10^-8

        # Sample the SD piece:
        lambda_omega_tpk = 1/sqrt(matrix(rgamma(n = (T-1)*(p-1)*K,
                                                shape = 1,
                                                rate = px_lambda_tpk + ss_omega_scale/2), nrow = T-1))
        px_lambda_tpk = matrix(rgamma(n = (T-1)*(p-1)*K,
                                      shape = 1,
                                      rate = 1/lambda_omega_tpk^2 + 1), nrow = T-1)
        #----------------------------------------------------------------------------
        # Second part: (jk)-specific shrinkage
        omega_scale = array(omega/lambda_omega_tpk, c(T-1, p-1, K))

        # Squared, then check for numerical issues:
        ss_omega_scale = colSums(omega_scale^2); ss_omega_scale = ss_omega_scale + (ss_omega_scale < 10^-16)*10^-8

        # Sample the SD piece:
        lambda_omega_pk = 1/sqrt(matrix(rgamma(n = (p-1)*K,
                                               shape = 1/2 + (T-1)/2,
                                               rate = px_lambda_pk + ss_omega_scale/2), nrow = p-1))
        px_lambda_pk = matrix(rgamma(n = (p-1)*K,
                                     shape = 1,
                                     rate = 1/lambda_omega_pk^2 + rep(1/lambda_omega_p^2, times = K)), nrow = p-1)
        #----------------------------------------------------------------------------
        # Third part: (j)-specific shrinkage
        lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                       shape = 1/2 + K/2,
                                       rate = px_lambda_p + rowSums(px_lambda_pk)))
        px_lambda_p = rgamma(n = p-1,
                             shape = 1,
                             rate = 1/lambda_omega_p^2 + 1/lambda_omega_0^2)
        #----------------------------------------------------------------------------
        # Fourth part: global shrinkage
        lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                       shape = 1/2 + (p-1)/2,
                                       rate = px_lambda_0 + sum(px_lambda_p)))
        px_lambda_0 = rgamma(n = 1,
                             shape = 1,
                             rate = 1/lambda_omega_0^2 + 1)
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters

        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep(matrix(lambda_omega_pk), each = (T-1))

        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      } else{

        stop("not yet fixed!")

        # Non-dynamic setting: grab the first one (all the same),
        # and store as (K x p-1) matrix
        #omega = t(alpha.arr[1, -1, ])
        omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

        #----------------------------------------------------------------------------
        # First part: jk-specific shrinkage

        # Scale by sigma_pk
        omega_scale = omega/t(sigma_omega_pk)

        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)

        # SD term:
        lambda_omega_kp = evolParams$sigma_wt

        #----------------------------------------------------------------------------
        # Second part: k-specific shrinkage (MGP)
        if(use_MGP){

          # Standardize, then reconstruct as matrix of size (T-1)*(p-1) x K:
          omega_scale = t(omega/lambda_omega_kp)

          # For each predictor:
          for(j in 1:(p-1)){
            delta_omega_pk[j,] = sampleMGP(theta.jh = omega_scale[j,],
                                           delta.h = delta_omega_pk[j,],
                                           a1 = a1_omega, a2 = a2_omega)
            tau_omega_pk[j,] = cumprod(delta_omega_pk[j,])
          }

          # Update the MGP parameters, common to all j:
          a1_omega = uni.slice(a1_omega, g = function(a){
            sum(dgamma(delta_omega_pk[,1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
          a2_omega = uni.slice(a2_omega, g = function(a){
            sum(dgamma(delta_omega_pk[,-1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
        }
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters

        # SD for omega:
        sigma_omega_pk = t(lambda_omega_kp)/sqrt(tau_omega_pk)

        # Cap at machine epsilon:
        sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(sigma_omega_pk^2))
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

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




# Modified dfosr_ar to include j-specific shrinkage:
dfosr_ar1 = function(Y, tau, X = NULL, K = NULL,
                    use_dynamic_reg = TRUE,
                    dist_reg_coef = "DHS", dist_reg_error = "NIG",
                    nsave = 1000, nburn = 1000, nskip = 10,
                    mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                    use_MGP = TRUE,
                    use_obs_SV = FALSE,
                    includeBasisInnovation = FALSE,
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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    nu = t(BtY) - BetaPsit; sigma_nu = sd(nu)
    theta = BetaPsit + nu; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_nu = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
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
  if(use_MGP || TRUE){
    a1_mu = 2; a2_mu = 3
    delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  } else {sigma_mu_k = abs(mu_k); px_mu_k = rep(1,K)}
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

    # MGP term:
    a1_omega = 2; a2_omega = 3
    delta_omega_k = rep(1,K); sigma_omega_k = 1/sqrt(cumprod(delta_omega_k))

    if(use_dynamic_reg){

      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])

      # Innovation:
      omega = diff(alpha_reg)

      # First, initialize the variable selection scaling terms:
      lambda_j = apply(array(omega, c(T-1, p-1, K)), 2, sd)
      lambda_0 = median(lambda_j) # Global term
      px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
      rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K) # Recurring term

      # Initialize the variance components:
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))
      # This is the scaling SD
      lambda_omega_tpk = evolParams$sigma_wt

      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j*rep(sigma_omega_k, each = (T-1)*(p-1))

      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      # Save the scaled omega for MGP computations below:
      if(use_MGP) omega_scale = matrix(array(omega/(lambda_omega_tpk*rep_lambda_j), c(T-1, p-1, K)), ncol = K)

    } else{

      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      #omega = t(alpha.arr[1, -1, ])
      omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

      # First, initialize the variable selection scaling terms:
      lambda_j = apply(omega, 2, sd)
      lambda_0 = median(lambda_j) # Global term
      px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
      rep_lambda_j = rep(lambda_j, each = K) # Recurring term

      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)

      # This is the scaling SD (Note: doesn't actually depend on T)
      lambda_omega_pk = evolParams$sigma_wt

      # SD for omega:
      sigma_omega_pk = lambda_omega_pk*rep(sigma_omega_k, times = (p-1))

      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))

      # Save the scaled omega for MGP computations below:
      if(use_MGP) omega_scale = t(omega/(lambda_omega_pk*rep_lambda_j))
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2 + sigma_nu^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    if(includeBasisInnovation){
      Y_tilde =  tcrossprod(theta, t(Psi)); sigma_tilde = sigma_nu
    } else {
      Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et
    }

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

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_nu^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_nu^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)

      sigma_nu = 1/sqrt(truncdist::rtrunc(1, "gamma",
                                          a = 10^-8, b = Inf,
                                          shape = (length(theta)+1)/2,
                                          rate = 1/2*sum((theta - BetaPsit)^2)))
    } else {theta = BetaPsit; sigma_nu = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}
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

    # Prior variance: MGP or Gaussian with half-Cauchy prior
    if(use_MGP || TRUE){
      # Mean Part
      delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                              delta.h = delta_mu_k,
                              a1 = a1_mu, a2 = a2_mu)
      sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

      # Variance part:
      # Standardize, then reconstruct as matrix of size T x K:
      delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                              delta.h = delta_eta_k,
                              a1 = a1_eta, a2 = a2_eta)
      sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
      a1_eta = uni.slice(a1_eta, g = function(a){
        dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
      a2_eta = uni.slice(a2_eta, g = function(a){
        sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

      # And standardize the residuals for future computations:
      eta_tk = eta_tk/rep(sigma_delta_k, each = T-1)
    } else {
      sigma_mu_k = 1/sqrt(rgamma(n = K, shape = 1/2 + 1/2, rate = mu_k^2/2 + px_mu_k))
      px_mu_k = rgamma(n = K, shape = 1/2 + 1/2, rate = 1/sigma_mu_k^2 + (1:K)^2)
    }

    # Sample the corresponding prior variance term(s):
    evolParams_int = sampleEvolParams(omega = eta_tk, evolParams =  evolParams_int, evol_error = dist_reg_error)
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
        # First part: j-specific shrinkage
        # Scale omega by sigma_k AND lambda_omega_tpk:
        omega_scale = omega/(lambda_omega_tpk*rep(sigma_omega_k, each = (T-1)*(p-1)))
        # Sum of squares (over t,k)
        ss_omega_j = apply(array(omega_scale^2, c(T-1, p-1, K)), 2, sum)
        # Offset:
        ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8

        # SD term:
        lambda_j = 1/sqrt(rgamma(n = p - 1,
                                 shape = 1/2 + (T-1)*K/2,
                                 rate = px_lambda_j + ss_omega_j/2))
        lambda_j = rep(1, p-1)
        # Sample the PX, global parameters outside this if-statement

        # Recurring term:
        rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K)
        #----------------------------------------------------------------------------
        # Second part: tjk-specific shrinkage

        # Scale by lambda_j AND sigma_k
        omega_scale = omega/(rep_lambda_j*rep(sigma_omega_k, each = (T-1)*(p-1)))

        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)

        # SD term:
        lambda_omega_tpk = evolParams$sigma_wt

        #----------------------------------------------------------------------------
        # Third part: k-specific shrinkage (MGP)
        if(use_MGP){
          omega_scale = matrix(array(omega/(lambda_omega_tpk*rep_lambda_j), c(T-1, p-1, K)), ncol = K)

          # Standardize, then reconstruct as matrix of size (T-1)*(p-1) x K:
          delta_omega_k = sampleMGP(theta.jh = omega_scale,
                                    delta.h = delta_omega_k,
                                    a1 = a1_omega, a2 = a2_omega)
          sigma_omega_k = 1/sqrt(cumprod(delta_omega_k))
          a1_omega = uni.slice(a1_omega, g = function(a){
            dgamma(delta_omega_k[1], shape = a, rate = 1, log = TRUE) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
          a2_omega = uni.slice(a2_omega, g = function(a){
            sum(dgamma(delta_omega_k[-1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)
        }
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters

        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j*rep(sigma_omega_k, each = (T-1)*(p-1))

        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      } else{

        # Non-dynamic setting: grab the first one (all the same),
        # and store as (K x p-1) matrix
        #omega = t(alpha.arr[1, -1, ])
        omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

        #----------------------------------------------------------------------------
        # First part: j-specific shrinkage
        # Scale omega by sigma_k AND lambda_omega_pk:
        omega_scale = omega/(lambda_omega_pk*rep(sigma_omega_k, times = (p-1)))
          # Sum of squares (over k)
        ss_omega_j = colSums(omega_scale^2)
          # Offset:
        ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8

        # SD term:
        lambda_j = 1/sqrt(rgamma(n = p - 1,
                                 shape = 1/2 + K/2,
                                 rate = px_lambda_j + ss_omega_j/2))
        lambda_j = rep(1, p-1)

        # Recurring term:
        rep_lambda_j = rep(lambda_j, each = K)
        #----------------------------------------------------------------------------
        # Second part: jk-specific shrinkage

        # Scale by lambda_j AND sigma_k
        omega_scale = omega/(rep_lambda_j*rep(sigma_omega_k, times = (p-1)))

        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)

        # SD term:
        lambda_omega_pk = evolParams$sigma_wt

        #----------------------------------------------------------------------------
        # Third part: k-specific shrinkage (MGP)
        if(use_MGP){

          # Standardize, then reconstruct as matrix of size (T-1)*(p-1) x K:
          omega_scale = t(omega/(lambda_omega_pk*rep_lambda_j))

          delta_omega_k = sampleMGP(theta.jh = omega_scale,
                                    delta.h = delta_omega_k,
                                    a1 = a1_omega, a2 = a2_omega)
          sigma_omega_k = 1/sqrt(cumprod(delta_omega_k))
          a1_omega = uni.slice(a1_omega, g = function(a){
            dgamma(delta_omega_k[1], shape = a, rate = 1, log = TRUE) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
          a2_omega = uni.slice(a2_omega, g = function(a){
            sum(dgamma(delta_omega_k[-1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)
        }
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters

        # SD for omega:
        sigma_omega_pk = lambda_omega_pk*rep_lambda_j*rep(sigma_omega_k, times = (p-1))

        # Cap at machine epsilon:
        sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))
      }
      #----------------------------------------------------------------------------
      # In both cases, need to sample the PX-parameters for lambda_j and global shrinkage parameters
      # These are the same for dynamic and non-dynamic cases
      # PX-j:
      px_lambda_j = rgamma(n = p - 1, shape = 1, rate = 1/lambda_j^2 + 1/lambda_0^2)

      # Global shrinkage:
      # SD term:
      lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_lambda_0 + sum(px_lambda_j)))

      # PX:
      px_lambda_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

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
# Modified dfosr_ar to include MGP priors for each predictor j:
dfosr_ar2 = function(Y, tau, X = NULL, K = NULL,
                     use_dynamic_reg = TRUE,
                     dist_reg_coef = "DHS", dist_reg_error = "NIG",
                     nsave = 1000, nburn = 1000, nskip = 10,
                     mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                     use_MGP = TRUE,
                     use_obs_SV = FALSE,
                     includeBasisInnovation = FALSE,
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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    nu = t(BtY) - BetaPsit; sigma_nu = sd(nu)
    theta = BetaPsit + nu; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_nu = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
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
  if(use_MGP){
    a1_mu = 2; a2_mu = 3
    delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  } else {sigma_mu_k = abs(mu_k); px_mu_k = rep(1,K)}
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

    # MGP term for each predictor:
    a1_omega = 2; a2_omega = 3
    delta_omega_pk = tau_omega_pk = matrix(1, nrow = p-1, ncol = K)

    if(use_dynamic_reg){

      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])

      # Innovation:
      omega = diff(alpha_reg)

      # Initialize the variance components:
      evolParams = initEvolParams(omega = omega, evol_error = dist_reg_coef)
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))

      # This is the scaling SD
      lambda_omega_tpk = evolParams$sigma_wt

      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep(matrix(1/sqrt(tau_omega_pk)), each = (T-1))

      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

    } else{

      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      #omega = t(alpha.arr[1, -1, ])
      omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

      evolParams = initEvolParams(omega = omega, evol_error = dist_reg_coef)

      # This is the scaling SD (Note: doesn't actually depend on T)
      lambda_omega_kp = evolParams$sigma_wt

      # SD for omega:
      sigma_omega_pk = t(lambda_omega_kp)/sqrt(tau_omega_pk)

      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(sigma_omega_pk^2))
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2 + sigma_nu^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    if(includeBasisInnovation){
      Y_tilde =  tcrossprod(theta, t(Psi)); sigma_tilde = sigma_nu
    } else {
      Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et
    }

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

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_nu^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_nu^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)

      sigma_nu = 1/sqrt(truncdist::rtrunc(1, "gamma",
                                          a = 10^-8, b = Inf,
                                          shape = (length(theta)+1)/2,
                                          rate = 1/2*sum((theta - BetaPsit)^2)))
    } else {theta = BetaPsit; sigma_nu = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}
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

    # Prior variance: MGP or Gaussian with half-Cauchy prior
    if(use_MGP){
      # Mean Part
      delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                              delta.h = delta_mu_k,
                              a1 = a1_mu, a2 = a2_mu)
      sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

      # Variance part:
      # Standardize, then reconstruct as matrix of size T x K:
      delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                              delta.h = delta_eta_k,
                              a1 = a1_eta, a2 = a2_eta)
      sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
      a1_eta = uni.slice(a1_eta, g = function(a){
        dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
      a2_eta = uni.slice(a2_eta, g = function(a){
        sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

      # And standardize the residuals for future computations:
      eta_tk = eta_tk/rep(sigma_delta_k, each = T-1)
    } else {
      sigma_mu_k = 1/sqrt(rgamma(n = K, shape = 1/2 + 1/2, rate = mu_k^2/2 + px_mu_k))
      px_mu_k = rgamma(n = K, shape = 1/2 + 1/2, rate = 1/sigma_mu_k^2 + (1:K)^2)
    }

    # Sample the corresponding prior variance term(s):
    evolParams_int = sampleEvolParams(omega = eta_tk, evolParams =  evolParams_int, evol_error = dist_reg_error)
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
        # First part: tjk-specific shrinkage

        # Scale by tau_omega_pk:
        omega_scale = omega*rep(matrix(sqrt(tau_omega_pk)), each = (T-1))

        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)

        # SD term:
        lambda_omega_tpk = evolParams$sigma_wt

        #----------------------------------------------------------------------------
        # Second part: k-specific shrinkage (MGP)
        if(use_MGP){
          # Scale the error term:
          omega_scale = array(omega/lambda_omega_tpk, c(T-1, p-1, K))

          # For each predictor:
          for(j in 1:(p-1)){
            delta_omega_pk[j,] = sampleMGP(theta.jh = omega_scale[,j,],
                                           delta.h = delta_omega_pk[j,],
                                           a1 = a1_omega, a2 = a2_omega)
            tau_omega_pk[j,] = cumprod(delta_omega_pk[j,])
          }

          # Update the MGP parameters, common to all j:
          a1_omega = uni.slice(a1_omega, g = function(a){
            sum(dgamma(delta_omega_pk[,1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
          a2_omega = uni.slice(a2_omega, g = function(a){
            sum(dgamma(delta_omega_pk[,-1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
        }
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters

        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep(matrix(1/sqrt(tau_omega_pk)), each = (T-1))

        # Cap at machine epsilon:
        sigma_omega_tpk[which(sigma_omega_tpk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      } else{

        # Non-dynamic setting: grab the first one (all the same),
        # and store as (K x p-1) matrix
        #omega = t(alpha.arr[1, -1, ])
        omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

        #----------------------------------------------------------------------------
        # First part: jk-specific shrinkage

        # Scale by sigma_pk
        omega_scale = omega/t(sigma_omega_pk)

        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)

        # SD term:
        lambda_omega_kp = evolParams$sigma_wt

        #----------------------------------------------------------------------------
        # Second part: k-specific shrinkage (MGP)
        if(use_MGP){

          # Standardize, then reconstruct as matrix of size (T-1)*(p-1) x K:
          omega_scale = t(omega/lambda_omega_kp)

          # For each predictor:
          for(j in 1:(p-1)){
            delta_omega_pk[j,] = sampleMGP(theta.jh = omega_scale[j,],
                                           delta.h = delta_omega_pk[j,],
                                           a1 = a1_omega, a2 = a2_omega)
            tau_omega_pk[j,] = cumprod(delta_omega_pk[j,])
          }

          # Update the MGP parameters, common to all j:
          a1_omega = uni.slice(a1_omega, g = function(a){
            sum(dgamma(delta_omega_pk[,1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
          a2_omega = uni.slice(a2_omega, g = function(a){
            sum(dgamma(delta_omega_pk[,-1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
        }
        #----------------------------------------------------------------------------
        # Last part:  update the variance parameters

        # SD for omega:
        sigma_omega_pk = t(lambda_omega_kp)/sqrt(tau_omega_pk)

        # Cap at machine epsilon:
        sigma_omega_pk[which(sigma_omega_pk < sqrt(.Machine$double.eps), arr.ind = TRUE)] = sqrt(.Machine$double.eps)

        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(sigma_omega_pk^2))
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

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


# Old code for AR:
dfosr_ar0 = function(Y, tau, X = NULL, K = NULL,
                    use_dynamic_reg = TRUE,
                    dist_reg_coef = "DHS", dist_reg_error = "NIG",
                    nsave = 1000, nburn = 1000, nskip = 10,
                    mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                    use_MGP = TRUE,
                    use_obs_SV = FALSE,
                    includeBasisInnovation = FALSE,
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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    nu = t(BtY) - BetaPsit; sigma_nu = sd(nu)
    theta = BetaPsit + nu; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_nu = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
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
  if(use_MGP){
    a1_mu = 2; a2_mu = 3
    delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  } else {sigma_mu_k = abs(mu_k); px_mu_k = rep(1,K)}
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

    # MGP term:
    a1_omega = 2; a2_omega = 3
    delta_omega_k = rep(1,K); sigma_omega_k = 1/sqrt(cumprod(delta_omega_k))

    if(use_dynamic_reg){

      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])

      evolParams = initEvolParams(omega = diff(alpha_reg), evol_error = dist_reg_coef)
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))

      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = evolParams$sigma_wt[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)
    } else{
      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      #omega = t(alpha.arr[1, -1, ])
      omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

      evolParams = initEvolParams(omega = omega, evol_error = dist_reg_coef)

      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(evolParams$sigma_wt^2)))
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2 + sigma_nu^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    if(includeBasisInnovation){
      Y_tilde =  tcrossprod(theta, t(Psi)); sigma_tilde = sigma_nu
    } else {
      Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et
    }

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

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_nu^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_nu^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)

      sigma_nu = 1/sqrt(truncdist::rtrunc(1, "gamma",
                                          a = 10^-8, b = Inf,
                                          shape = (length(theta)+1)/2,
                                          rate = 1/2*sum((theta - BetaPsit)^2)))
    } else {theta = BetaPsit; sigma_nu = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}
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

    # Prior variance: MGP or Gaussian with half-Cauchy prior
    if(use_MGP){
      # Mean Part
      delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                              delta.h = delta_mu_k,
                              a1 = a1_mu, a2 = a2_mu)
      sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

      # Variance part:
      # Standardize, then reconstruct as matrix of size T x K:
      delta_eta_k = sampleMGP(theta.jh = matrix(eta_tk/evolParams_int$sigma_wt, ncol = K),
                              delta.h = delta_eta_k,
                              a1 = a1_eta, a2 = a2_eta)
      sigma_delta_k = 1/sqrt(cumprod(delta_eta_k))
      a1_eta = uni.slice(a1_eta, g = function(a){
        dgamma(delta_eta_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
      a2_eta = uni.slice(a2_eta, g = function(a){
        sum(dgamma(delta_eta_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

      # And standardize the residuals for future computations:
      eta_tk = eta_tk/rep(sigma_delta_k, each = T-1)
    } else {
      sigma_mu_k = 1/sqrt(rgamma(n = K, shape = 1/2 + 1/2, rate = mu_k^2/2 + px_mu_k))
      px_mu_k = rgamma(n = K, shape = 1/2 + 1/2, rate = 1/sigma_mu_k^2 + (1:K)^2)
    }

    # Sample the corresponding prior variance term(s):
    evolParams_int = sampleEvolParams(omega = eta_tk, evolParams =  evolParams_int, evol_error = dist_reg_error)
    evolParams0_int = sampleEvol0(mu0 = as.matrix(gamma_tk[1,]), evolParams0 = evolParams0_int)

    # Update the error SD for gamma:
    sigma_eta_tk = evolParams_int$sigma_wt*rep(sigma_delta_k, each = T-1)

    # Evolution and initial variance:
    for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_eta_tk[,k]^2
    diag(W0[ind.intercept, ind.intercept]) = as.numeric(evolParams0_int$sigma_w0^2)
    #----------------------------------------------------------------------------
    # Step 7: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){

      if(use_dynamic_reg){ # Dynamic setting

        alpha_reg = as.matrix(alpha[,-ind.intercept])

        # Innovations:
        omega = diff(alpha_reg)

        # MGP, if desirable:
        if(use_MGP){
          # Standardize, then reconstruct as matrix of size (T-1)*(p-1) x K:
          delta_omega_k = sampleMGP(theta.jh = matrix(array(omega/evolParams$sigma_wt, c(T-1, p-1, K)), ncol = K),
                                    delta.h = delta_omega_k,
                                    a1 = a1_omega, a2 = a2_omega)
          sigma_omega_k = 1/sqrt(cumprod(delta_omega_k))
          a1_omega = uni.slice(a1_omega, g = function(a){
            dgamma(delta_omega_k[1], shape = a, rate = 1, log = TRUE) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
          a2_omega = uni.slice(a2_omega, g = function(a){
            sum(dgamma(delta_omega_k[-1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

          # And standardize the residuals:
          omega = omega/rep(sigma_omega_k, each = (T-1)*(p-1))
        }

        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)

        # Update the variances:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = evolParams$sigma_wt[,j]^2*(rep(sigma_omega_k^2, each = p-1)[j])
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      } else{
        # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
        #omega = t(alpha.arr[1, -1, ])
        omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

        # MGP, if desirable:
        if(use_MGP){
          # Standardize, then reconstruct as matrix of size (p-1) x K:
          delta_omega_k = sampleMGP(theta.jh = t(omega/evolParams$sigma_wt),
                                    delta.h = delta_omega_k,
                                    a1 = a1_omega, a2 = a2_omega)
          sigma_omega_k = 1/sqrt(cumprod(delta_omega_k))
          a1_omega = uni.slice(a1_omega, g = function(a){
            dgamma(delta_omega_k[1], shape = a, rate = 1, log = TRUE) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
          a2_omega = uni.slice(a2_omega, g = function(a){
            sum(dgamma(delta_omega_k[-1], shape = a, rate = 1, log = TRUE)) +
              dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)

          # And standardize the residuals:
          omega = omega/rep(sigma_omega_k, times = p-1)
        }

        # And sample the parameters:
        evolParams = sampleEvolParams(omega = omega, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)

        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(evolParams$sigma_wt^2)))*(rep(sigma_omega_k^2, each = p-1))
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

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

# Independent case
dfosr_ind0 = function(Y, tau, X = NULL, K = NULL,
                     use_dynamic_reg = TRUE,
                     dist_reg_coef = "DHS", dist_reg_error = "NIG",
                     nsave = 1000, nburn = 1000, nskip = 10,
                     mcmc_params = list("beta", "fk", "alpha", "mu_k", "ar_phi"),
                     use_MGP = TRUE,
                     use_obs_SV = FALSE,
                     includeBasisInnovation = FALSE,
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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    nu = t(BtY) - BetaPsit; sigma_nu = sd(nu)
    theta = BetaPsit + nu; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_nu = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
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
  # MGP SD: multiplicative, so fix at 1 (in case we don't use it)
  sigma_k = rep(1, K) # SD term
  if(use_MGP){
    a1 = 2; a2 = 3 # Hyperparameters
    delta_k = rep(1, K) # Product term
  }

  # Overall mean term (and T x K case)
  mu_k = as.matrix(colMeans(Beta)); mu_tk = matrix(rep(mu_k, each =  T), nrow = T)

  # SD for mu_k:
  lambda_mu = sd(mu_k); px_mu = 1
  sigma_mu_k = sigma_k*lambda_mu  # This might be redundant throughout...
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

  # Intercept (or gamma) components:
  gamma_tk = as.matrix(alpha[,ind.intercept])

  # Initialize the corresponding SD term(s):
  evolParams_int = initEvolParams(omega = gamma_tk, evol_error = dist_reg_error)
  # This is the scaling SD
  lambda_gamma_tk = evolParams_int$sigma_wt
  # SD for gamma_tk:
  sigma_gamma_tk = lambda_gamma_tk*rep(sigma_k, each = T)

  # Evolution and initial variance:
  for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_gamma_tk[-1, k]^2
  diag(W0[ind.intercept, ind.intercept]) = as.numeric(sigma_gamma_tk[1,]^2)
  #----------------------------------------------------------------------------
  # Non-intercept term:
  if(p > 1){

    if(use_dynamic_reg){
      # Dynamic setting:
      alpha_reg = as.matrix(alpha[,-ind.intercept])

      # Innovation:
      omega = diff(alpha_reg)

      # First, initialize the variable selection scaling terms:
      lambda_j = apply(array(omega, c(T-1, p-1, K)), 2, sd)
      lambda_0 = median(lambda_j) # Global term
      px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
      rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K) # Recurring term

      # Initialize the variance components:
      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)
      evolParams0 = initEvol0(mu0 = as.matrix(alpha_reg[1,]))
      # This is the scaling SD
      lambda_omega_tpk = evolParams$sigma_wt

      # SD for omega:
      sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j*rep(sigma_k, each = (T-1)*(p-1))

      # Evolution and initial variance:
      for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      # Save the scaled omega for MGP computations below:
      if(use_MGP) omega_scale = matrix(array(omega/(lambda_omega_tpk*rep_lambda_j), c(T-1, p-1, K)), ncol = K)

    } else{
      # Non-dynamic setting: grab the first one (all the same) and store as (K x p-1) matrix
      #omega = t(alpha.arr[1, -1, ])
      omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

      # First, initialize the variable selection scaling terms:
      lambda_j = apply(omega, 2, sd)
      lambda_0 = median(lambda_j) # Global term
      px_lambda_j = rep(1, p-1); px_lambda_0 = 1 # PX terms
      rep_lambda_j = rep(lambda_j, each = K) # Recurring term

      evolParams = initEvolParams(omega = omega/rep_lambda_j, evol_error = dist_reg_coef)

      # This is the scaling SD (Note: doesn't actually depend on T)
      lambda_omega_pk = evolParams$sigma_wt

      # SD for omega:
      sigma_omega_pk = lambda_omega_pk*rep(sigma_k, times = (p-1))

      # Only need the inital variance in the static case:
      diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))

      # Save the scaled omega for MGP computations below:
      if(use_MGP) omega_scale = t(omega/(lambda_omega_pk*rep_lambda_j))
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2 + sigma_nu^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance depend on basis innovation:
    if(includeBasisInnovation){
      Y_tilde =  tcrossprod(theta, t(Psi)); sigma_tilde = sigma_nu
    } else {
      Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et
    }

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

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)
    #----------------------------------------------------------------------------
    # Step 4: Sample the basis terms (if desired)
    #----------------------------------------------------------------------------
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_nu^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_nu^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_nu^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)

      sigma_nu = 1/sqrt(truncdist::rtrunc(1, "gamma",
                                          a = 10^-8, b = Inf,
                                          shape = (length(theta)+1)/2,
                                          rate = 1/2*sum((theta - BetaPsit)^2)))
    } else {theta = BetaPsit; sigma_nu = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}
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
    postSD = 1/sqrt(colSums(1/sigma_gamma_tk^2) + 1/sigma_mu_k^2)
    postMean = (colSums(gamma_tk_c/sigma_gamma_tk^2))*postSD^2
    mu_k = rnorm(n = K, mean = postMean, sd = postSD)
    mu_tk = matrix(rep(mu_k, each =  T), nrow = T)

    # And update the non-centered parameter:
    gamma_tk = gamma_tk_c - mu_tk
    #----------------------------------------------------------------------------
    # Step 7: Variance Components (non-MGP) + aux variables
    #----------------------------------------------------------------------------
    # 7a: Variance Components for mu_k
    # SD term:
    lambda_mu = 1/sqrt(rgamma(n = 1, shape = K/2 + 1/2, rate = sum((mu_k/sigma_k)^2)/2 + px_mu))
    px_mu = rgamma(n = 1, shape = 1, rate = 1/lambda_mu^2 + 1/100^2)
    #----------------------------------------------------------------------------
    # 7b: Variance Components for gamma_tk
    evolParams_int = sampleEvolParams(omega = gamma_tk/rep(sigma_k, each = T),
                                      evolParams =  evolParams_int,
                                      evol_error = dist_reg_error)
    # SD term:
    lambda_gamma_tk = evolParams_int$sigma_wt
    #----------------------------------------------------------------------------
    # 7c Variance components for omega_jkt
    if(p > 1){

      if(use_dynamic_reg){

        # Dynamic setting

        # Regression (non-intercept) coefficients
        alpha_reg = as.matrix(alpha[,-ind.intercept])

        # Random walk, so compute difference for innovations:
        omega = diff(alpha_reg)

        #----------------------------------------------------------------------------
        # First part: j-specific shrinkage
        # Scale omega by sigma_k AND lambda_omega_tpk:
        omega_scale = omega/(lambda_omega_tpk*rep(sigma_k, each = (T-1)*(p-1)))
        # Sum of squares (over t,k)
        ss_omega_j = apply(array(omega_scale^2, c(T-1, p-1, K)), 2, sum)
        # Offset:
        ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8

        # SD term:
        lambda_j = 1/sqrt(rgamma(n = p - 1,
                                 shape = 1/2 + (T-1)*K/2,
                                 rate = px_lambda_j + ss_omega_j/2))
        #lambda_j = rep(1, p-1)
        # Sample the PX, global parameters outside this if-statement

        # Recurring term:
        rep_lambda_j = rep(rep(lambda_j, each = T-1), times = K)
        #----------------------------------------------------------------------------
        # Second part: tjk-specific shrinkage

        # Scale by lambda_j AND sigma_k
        omega_scale = omega/(rep_lambda_j*rep(sigma_k, each = (T-1)*(p-1)))

        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)
        evolParams0 = sampleEvol0(mu0 = as.matrix(alpha_reg[1,]), evolParams0, A = 1)

        # SD term:
        lambda_omega_tpk = evolParams$sigma_wt

        # Save the scaled omega for MGP computations later:
        if(use_MGP) omega_scale = matrix(array(omega/(lambda_omega_tpk*rep_lambda_j), c(T-1, p-1, K)), ncol = K)

      } else {

        # Non-dynamic setting: grab the first one (all the same),
        # and store as (K x p-1) matrix
        #omega = t(alpha.arr[1, -1, ])
        omega = matrix(t(alpha.arr[1, -1, ]), nrow = K)

        #----------------------------------------------------------------------------
        # First part: j-specific shrinkage
        # Scale omega by sigma_k AND lambda_omega_pk:
        omega_scale = omega/(lambda_omega_pk*rep(sigma_k, times = (p-1)))
        # Sum of squares (over k)
        ss_omega_j = colSums(omega_scale^2)
        # Offset:
        ss_omega_j = ss_omega_j + (ss_omega_j < 10^-16)*10^-8

        # SD term:
        lambda_j = 1/sqrt(rgamma(n = p - 1,
                                 shape = 1/2 + K/2,
                                 rate = px_lambda_j + ss_omega_j/2))
        #lambda_j = rep(1, p-1)

        # Recurring term:
        rep_lambda_j = rep(lambda_j, each = K)
        #----------------------------------------------------------------------------
        # Second part: jk-specific shrinkage

        # Scale by lambda_j AND sigma_k
        omega_scale = omega/(rep_lambda_j*rep(sigma_k, times = (p-1)))

        # Sample the parameters
        evolParams = sampleEvolParams(omega = omega_scale, evolParams, 1/sqrt(T), evol_error = dist_reg_coef)

        # SD term:
        lambda_omega_pk = evolParams$sigma_wt

        # Save the scaled omega for MGP computations later:
        if(use_MGP) omega_scale = t(omega/(lambda_omega_pk*rep_lambda_j))
      }
      #----------------------------------------------------------------------------
      # Lastly, sample the PX-parameters for lambda_j and global shrinkage parameters
      # These are the same for dynamic and non-dynamic cases
      # PX-j:
      px_lambda_j = rgamma(n = p - 1, shape = 1, rate = 1/lambda_j^2 + 1/lambda_0^2)

      # Global shrinkage:
      # SD term:
      lambda_0 = 1/sqrt(rgamma(n = 1, shape = 1/2 + (p-1)/2, rate = px_lambda_0 + sum(px_lambda_j)))

      # PX:
      px_lambda_0 = rgamma(n = 1, shape = 1, rate = 1/lambda_0^2 + 1)
    }

    #----------------------------------------------------------------------------
    # Step 8: Sample the MGP parameters
    #----------------------------------------------------------------------------
    if(use_MGP){
      if(p > 1){
        vecMGP = rbind(
          matrix(mu_k/lambda_mu, ncol = K),
          matrix(gamma_tk/lambda_gamma_tk, ncol = K),
          matrix(omega_scale, ncol = K))
      } else vecMGP = rbind(
        matrix(mu_k/lambda_mu, ncol = K),
        matrix(gamma_tk/lambda_gamma_tk, ncol = K))

      # Product terms:
      delta_k =  sampleMGP(theta.jh = vecMGP,delta.h = delta_k,a1 = a1, a2 = a2)

      # Standard deviation terms:
      sigma_k = 1/sqrt(cumprod(delta_k))

      # Option to sample the hyperparameters, a1 and a2:
      a1 = uni.slice(a1, g = function(a){
        dgamma(delta_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = 100)
      a2 = uni.slice(a2, g = function(a){
        sum(dgamma(delta_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = 100)
    }
    #----------------------------------------------------------------------------
    # Step 9: Update the TOTAL variance terms
    #----------------------------------------------------------------------------

    # Variance of the mean term(s):
    sigma_mu_k = sigma_k*lambda_mu

    # Variance of the factor innovation (intercept)
    sigma_gamma_tk = lambda_gamma_tk*rep(sigma_k, each = T)

    # For updating the KFAS object later:
    for(k in 1:K) Wt[ind.intercept, ind.intercept,][k,k,-T] = sigma_gamma_tk[-1, k]^2
    diag(W0[ind.intercept, ind.intercept]) = as.numeric(sigma_gamma_tk[1,]^2)

    # And the regression terms:
    if(p > 1){
      if(use_dynamic_reg){
        # SD for omega:
        sigma_omega_tpk = lambda_omega_tpk*rep_lambda_j*rep(sigma_k, each = (T-1)*(p-1))

        # Evolution and initial variance:
        for(j in 1:(K*(p-1))) Wt[-ind.intercept, -ind.intercept,][j,j,-T] = sigma_omega_tpk[,j]^2
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(evolParams0$sigma_w0^2)

      } else{
        # SD for omega:
        sigma_omega_pk = lambda_omega_pk*rep_lambda_j*rep(sigma_k, times = (p-1))

        # Only need the inital variance in the static case:
        diag(W0[-ind.intercept, -ind.intercept]) = as.numeric(matrix(t(sigma_omega_pk^2)))
      }
    }
    #----------------------------------------------------------------------------
    # Step 10: Adjust the ordering
    #----------------------------------------------------------------------------
    if(nsi == 10 && K > 1){adjOrder = order(tau_f_k, decreasing = TRUE); tau_f_k = tau_f_k[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

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
