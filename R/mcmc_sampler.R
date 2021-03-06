#' @useDynLib fosr
#' @importFrom Rcpp sourceCpp
NULL

#' MCMC Sampling Algorithm for the Function-on-Scalars Regression Model
#'
#' Runs the MCMC for the function-on-scalars regression model based on
#' an FDLM-type expansion. Here we assume the factor regression has independent errors,
#' which allows for subject-specific random effects,
#' as well as some additional default conditions.
#'
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of subjects and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{n x p} matrix of predictors; if NULL, only include an intercept
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
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
#' \item "sigma_e" (observation error SD)
#' \item "sigma_g" (random effects SD)
#' \item "Yhat" (fitted values)
#' \item "trsigma" (the sum of trace Sigma_i)
#' \item "sigma_delta_k" (random effects gamma specific SD)
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
#' out = fosr(Y = Y, tau = tau, X = X, K = 6, mcmc_params = list("fk", "alpha", "Yhat"))
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
#' @import truncdist
#' @export
fosr = function(Y, tau, X = NULL, K = NULL,
                nsave = 1000, nburn = 1000, nskip = 3,
                mcmc_params = list("beta", "fk", "alpha", "sigma_e", "sigma_g", "trsigma"),
                computeDIC = TRUE){

  # Some options (for now):
  sample_nu = TRUE # Sample DF parameter, or fix at nu=3?
  sample_a1a2 = TRUE # Sample a1, a2, or fix at a1=2, a2=3?

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

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init_d(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  K = ncol(Beta) # to be sure we have the right value

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing){na.ind = which(is.na(Yna), arr.ind = TRUE); Y = inits$Y0}
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  Yhat = tcrossprod(Beta, Fmat)

  # Initialize the (time-dependent) observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)

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
  a1_mu = 2; a2_mu = 3
  delta_mu_k = sampleMGP(matrix(mu_k, ncol = K), rep(1,K), a1 = a1_mu, a2 = a2_mu)
  sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
  #----------------------------------------------------------------------------
  # Initialize the corresponding SD term(s):
  xi_gamma_ik = 1/gamma_ik^2; # Precision scale
  nu = 3  # (initial) degrees of freedom

  # MGP term:
  a1_gamma = 2; a2_gamma = 3;
  delta_gamma_k = rep(1,K); sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))

  # Update the error SD for gamma:
  sigma_gamma_ik = rep(sigma_delta_k, each = n)/sqrt(xi_gamma_ik)
  #----------------------------------------------------------------------------
  if(p > 1){
    omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

    # predictor p, factor k:
    sigma_omega_pk = abs(omega)
    xi_omega_pk = matrix(1, nrow = p-1, ncol = K) # PX term

    # predictor p:
    lambda_omega_p = rowMeans(sigma_omega_pk)
    xi_omega_p = rep(1, (p-1)) # PX term

    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
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
  if(!is.na(match('trsigma', mcmc_params))) post.trsigma = array(NA, c(nsave))
  if(!is.na(match('sigma_delta_k', mcmc_params))) post.sigma_delta_k = array(NA, c(nsave, K))
  if(!is.na(match('nu', mcmc_params))) post.nu = array(NA, c(nsave))

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
    if(any.missing){
      Y[na.ind] = Yhat[na.ind] + sigma_et[na.ind[,1]]*rnorm(nrow(na.ind))
      BtY = tcrossprod(t(splineInfo$Bmat), Y)
    }
    #----------------------------------------------------------------------------
    # Step 2: Sample the FLCs
    #----------------------------------------------------------------------------
    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = splineInfo$BtB, #diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;

    # Sample the smoothing parameters:
    #tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = TRUE)
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = FALSE)
    #----------------------------------------------------------------------------
    # Step 3: Sample the regression coefficients (and therefore the factors)
    #----------------------------------------------------------------------------
    # Pseudo-response and pseudo-variance:
    Y_tilde = crossprod(BtY, Psi); sigma_tilde = sigma_et

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

    # Prior variance: MGP
    # Mean Part
    delta_mu_k =  sampleMGP(theta.jh = matrix(mu_k, ncol = K),
                            delta.h = delta_mu_k,
                            a1 = a1_mu, a2 = a2_mu)
    sigma_mu_k = 1/sqrt(cumprod(delta_mu_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_mu = uni.slice(a1_mu, g = function(a){
        dgamma(delta_mu_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_mu = uni.slice(a2_mu,g = function(a){
        sum(dgamma(delta_mu_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Variance part:
    # Standardize, then reconstruct as matrix of size n x K:
    delta_gamma_k = sampleMGP(theta.jh = matrix(gamma_ik*sqrt(xi_gamma_ik), ncol = K),
                              delta.h = delta_gamma_k,
                              a1 = a1_gamma, a2 = a2_gamma)
    sigma_delta_k = 1/sqrt(cumprod(delta_gamma_k))
    # And hyperparameters:
    if(sample_a1a2){
      a1_gamma = uni.slice(a1_gamma, g = function(a){
        dgamma(delta_gamma_k[1], shape = a, rate = 1, log = TRUE) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)}, lower = 0, upper = Inf)
      a2_gamma = uni.slice(a2_gamma, g = function(a){
        sum(dgamma(delta_gamma_k[-1], shape = a, rate = 1, log = TRUE)) +
          dgamma(a, shape = 2, rate = 1, log = TRUE)},lower = 0, upper = Inf)
    }

    # Sample the corresponding prior variance term(s):
    xi_gamma_ik = matrix(rgamma(n = n*K,
                                shape = nu/2 + 1/2,
                                rate = nu/2 + (gamma_ik/rep(sigma_delta_k, each = n))^2/2), nrow = n)
    # Sample degrees of freedom?
    if(sample_nu){
      nu = uni.slice(nu, g = function(nu){
        sum(dgamma(xi_gamma_ik, shape = nu/2, rate = nu/2, log = TRUE)) +
          dunif(nu, min = 2, max = 128, log = TRUE)}, lower = 2, upper = 128)
    }

    # Update the error SD for gamma:
    sigma_gamma_ik = rep(sigma_delta_k, each = n)/sqrt(xi_gamma_ik)
    #----------------------------------------------------------------------------
    # Step 6: Sample the non-intercept parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){

      omega = matrix(alpha_pk[-1,], nrow = p-1) # Not the intercept

      #----------------------------------------------------------------------------
      # predictor p, factor k:
      omega2 = omega^2; omega2 = omega2 + (omega2 < 10^-16)*10^-8
      sigma_omega_pk = matrix(1/sqrt(rgamma(n = (p-1)*K,
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_pk + omega2/2)), nrow = p-1)
      xi_omega_pk = matrix(rgamma(n = (p-1)*K,
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, times = K) + 1/sigma_omega_pk^2), nrow = p-1)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + K/2,
                                     rate = xi_omega_p + rowSums(xi_omega_pk)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Step 7: Adjust the ordering
    #----------------------------------------------------------------------------
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
        if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('sigma_g', mcmc_params))) post.sigma_g[isave,,] = sigma_gamma_ik
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Yhat # + sigma_e*rnorm(length(Y))
        if(!is.na(match('trsigma', mcmc_params))) post.trsigma[isave] = n*m*sigma_e^2 + sum(sigma_gamma_ik^2)
        if(!is.na(match('sigma_delta_k', mcmc_params))) post.sigma_delta_k[isave,] = sigma_delta_k
        if(!is.na(match('nu', mcmc_params))) post.nu[isave] = nu
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
  if(!is.na(match('trsigma', mcmc_params))) mcmc_output$trsigma = post.trsigma
  if(!is.na(match('sigma_delta_k', mcmc_params))) mcmc_output$sigma_delta_k = post.sigma_delta_k
  if(!is.na(match('nu', mcmc_params))) mcmc_output$nu = post.nu

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
