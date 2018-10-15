source("../R/helper_functions.R")

#' MCMC Sampling Algorithm for the Function-on-Scalars Regression Model
#' with pre-fixed basis functions
#'
#' Runs the MCMC for the function-on-scalars regression model based on
#' a known basis expansion. Here we assume the factor regression has independent errors,
#' which allows for subject-specific random effects,
#' as well as some additional default conditions.
#'
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of subjects and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{n x p} matrix of predictors; if NULL, only include an intercept
#' @param use_fpca logical; if TRUE, use FPCA, otherwise use orthonormal B-splines
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
#' \item "sigma_e" (observation error SD)
#' \item "sigma_g" (random effects SD)
#' \item "Yhat" (fitted values)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{nm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x n x M},  may be inefficient
#'
#' @examples
#' # Simulate some data:
#' sim_data = simulate_fosr(n = 100, m = 20, p_0 = 100, p_1 = 5)
#'
#' # Data:
#' Y = sim_data$Y; X = sim_data$X; tau = sim_data$tau
#'
#' # Dimensions:
#' n = nrow(Y); m = ncol(Y); p = ncol(X)
#'
#' # Run the FOSR:
#' out = fosr_basis(Y = Y, tau = tau, X = X, mcmc_params = list("fk", "alpha", "Yhat"))
#'
#' # Plot a posterior summary of a regression function, say j = 3:
#' j = 3; post_alpha_tilde_j = get_post_alpha_tilde(out$fk, out$alpha[,j,])
#' plot_curve(post_alpha_tilde_j, tau = tau)
#' # Add the true curve:
#' lines(tau, sim_data$alpha_tilde_true[,j], lwd=6, col='green', lty=6)
#'
#' # Plot the loading curves:
#' plot_flc(out$fk, tau = tau)
#'
#' # Plot the fitted values for a random subject:
#' i = sample(1:n, 1)
#' plot_fitted(y = Y[i,], mu = colMeans(out$Yhat[,i,]),
#'             postY = out$Yhat[,i,], y_true = sim_data$Y_true[i,], t01 = tau)
#'
#' @import truncdist refund
#' @export
fosr_basis = function(Y, tau, X = NULL, use_fpca = FALSE,
                nsave = 1000, nburn = 1000, nskip = 3,
                mcmc_params = list("beta", "fk", "alpha", "sigma_e", "sigma_g"),
                computeDIC = TRUE){
  #----------------------------------------------------------------------------
  # Assume that we've done checks elsewhere
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  n = nrow(Y); m = ncol(Y); d = ncol(tau)

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

  # Initialize the factors
  Beta = tcrossprod(Y, t(Fmat))
  K = ncol(Beta) # to be sure we have the right value

  # Initialize the conditional expectation:
  Yhat = tcrossprod(Beta, Fmat)

  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)
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
  X = cbind(rep(1, n), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')

  # Number of predictors (including the intercept)
  p = ncol(X)
  #----------------------------------------------------------------------------
  # Initialize the regression terms (and the mean term)
  alpha_pk = matrix(0, nrow = p, ncol = K) # Regression coefficients
  gamma_ik = matrix(0, nrow = n, ncol = K) # Residuals

  # Initialize the regression coefficients via sampling (p >= n) or OLS (p < n)
  for(k in 1:K) {
    if(p >= n){
      alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_et,
                                        Ddiag = rep(.01*sigma_e^2, p),
                                        alpha = tcrossprod(Y, t(Fmat[,k]))/sigma_e)
    } else alpha_pk[,k] = lm(Beta[,k] ~ X - 1)$coef

    # Residuals:
    gamma_ik[,k] = Beta[,k] - X%*%alpha_pk[,k]
  }

  # Intercept term:
  mu_k = as.matrix(alpha_pk[1,])

  # SD term for mu_k:
  sigma_mu_k = rep(100, K) # This will stay fixed
  #----------------------------------------------------------------------------
  # SD term for gamma_ik:
  sigma_gamma_ik = matrix(rep(apply(gamma_ik, 2, sd), each = n), nrow = n)
  #----------------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept
    sigma_omega_pk= matrix(rep(apply(omega, 1, sd), times = K), nrow = p-1)
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, n, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, p, K))
  if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_g', mcmc_params))) post.sigma_g = array(NA, c(nsave, n, K))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, n, m))
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

    # Not necessary here

    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------

    # Not necessary here

    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = tcrossprod(Y, t(Fmat)); sigma_tilde = sigma_et

    # Draw Separately for each k:
    for(k in 1:K){
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k]; sigma_tilde_k = sqrt(sigma_tilde^2 + sigma_gamma_ik[,k]^2)

      if(p >= n){
        # Fast sampler for p >= n (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(Phi = X/sigma_tilde_k,
                                          Ddiag = as.numeric(c(sigma_mu_k[k],sigma_omega_pk[,k])^2),
                                          alpha = y_tilde_k/sigma_tilde_k)
      } else {
        # Fast sampler for p < n (Rue, 2001?)
        if(p > 1){
          chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k],sigma_omega_pk[,k])^2)))
        } else chQ_k = chol(crossprod(X/sigma_tilde_k) + diag(as.numeric(1/c(sigma_mu_k[k])^2), p))
        ell_k = crossprod(X, y_tilde_k/sigma_tilde_k^2)
        alpha_pk[,k] = backsolve(chQ_k, forwardsolve(t(chQ_k), ell_k) + rnorm(p))
      }
    }

    # And sample the errors gamma_ik:
    postSD = 1/sqrt(rep(1/sigma_tilde^2, times = K) + matrix(1/sigma_gamma_ik^2))
    postMean = matrix((Y_tilde - X%*%alpha_pk)/rep(sigma_tilde^2, times = K))*postSD^2
    gamma_ik = matrix(rnorm(n = n*K, mean = postMean, sd = postSD), nrow = n)

    # Update the factors:
    Beta = X%*%alpha_pk + gamma_ik

    # And the fitted curves:
    Yhat = tcrossprod(Beta, Fmat)
    #----------------------------------------------------------------------------
    # Step 4: Sample the observation error variance
    #----------------------------------------------------------------------------
    # Or use uniform prior?
    sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, n)

    #----------------------------------------------------------------------------
    # Step 5: Sample the intercept/gamma parameters (Note: could use ASIS)
    #----------------------------------------------------------------------------
    mu_k = alpha_pk[1,]

    # SD term for gamma_ik:
    sigma_gamma_ik = matrix(rep(apply(gamma_ik, 2, function(x){
      1/sqrt(rgamma(n = 1, shape = 0.001 + n/2, rate = 0.001 + sum(x^2)/2))
    }), each = n), nrow = n)
    #----------------------------------------------------------------------------

    #----------------------------------------------------------------------------
    # Step 6: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      # SD term for omega_pk:
      sigma_omega_pk = matrix(rep(apply(omega, 1, function(x){
        1/sqrt(rgamma(n = 1, shape = 0.001 + K/2, rate = 0.001 + sum(x^2)/2))
      }), times = K), nrow = p-1)
    }
    #----------------------------------------------------------------------------
    # Step 7: Adjust the ordering
    #----------------------------------------------------------------------------
    # Not necessary here

    #----------------------------------------------------------------------------
    # Store the MCMC output:
    #----------------------------------------------------------------------------
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
        if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('sigma_g', mcmc_params))) post.sigma_g[isave,,] = sigma_gamma_ik
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Yhat # + sigma_e*rnorm(length(Y))
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Yhat), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('sigma_g', mcmc_params))) mcmc_output$sigma_g = post.sigma_g
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_e), m*n),
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
