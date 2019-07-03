#----------------------------------------------------------------------------
#' Simulate a function-on-scalar regression model
#'
#' Simulate data from a function-on-scalar regression model, allowing for
#' subject-specific random effects. The predictors are multivariate normal with
#' mean zero and covariance \code{corr^abs(j1-j2)} for correlation parameter \code{corr}
#' between predictors \code{j1} and \code{j2}.
#' More predictors than observations (p > n) is allowed.
#'
#' @param n number of observed curves (i.e., number of subjects)
#' @param m total number of observation points (i.e., points along the curve)
#' @param RSNR root signal-to-noise ratio
#' @param K_true rank of the model (i.e., number of basis functions used for the functional data simulations)
#' @param p_0 number of true zero regression coefficients
#' @param p_1 number of true nonzero regression coefficients
#' @param sparse_factors logical; if TRUE, then for each nonzero predictor j,
#' sample a subset of k=1:K_true factors to be nonzero#'
#' @param corr correlation parameter for predictors
#' @param perc_missing percentage of missing data (between 0 and 1); default is zero
#' @param X the design matrix to use. By default, \code{X} is NULL. When \code{X} is NULL, a design matrix
#'          is created. If \code{X} is specified, then params \code{p_0}, \code{n} and
#           \code{corr} are not used.
#'
#' @return a list containing the following:
#' \itemize{
#' \item \code{Y}: the simulated \code{n x m} functional data matrix
#' \item \code{X}: the simulated \code{n x p} design matrix
#' \item \code{tau}: the \code{m}-dimensional vector of observation points
#' \item \code{Y_true}: the true \code{n x m} functional data matrix (w/o noise)
#' \item \code{alpha_tilde_true} the true \code{m x p} matrix of regression coefficient functions
#' \item \code{alpha_arr_true} the true \code{K_true x p} matrix of regression coefficient factors
#' \item \code{Beta_true} the true \code{n x K_true} matrix of factors
#' \item \code{F_true} the true \code{m x K_true} matrix of basis (loading curve) functions
#' \item \code{sigma_true} the true observation error standard deviation
#' }
#'
#' @note The basis functions (or loading curves) are orthonormalized polynomials,
#' so large values of \code{K_true} are not recommended.
#'
#' @examples
#' # Example: simulate FOSR
#' sim_data = simulate_fosr(n = 100, m = 20, p_0 = 100, p_1 = 5)
#' Y = sim_data$Y; X = sim_data$X; tau = sim_data$tau
#'
#' @import truncdist
#' @export
simulate_fosr = function(n = 100,
                         m = 50,
                         RSNR = 5,
                         K_true = 4,
                         p_0 = 1000,
                         p_1 = 5,
                         sparse_factors = TRUE,
                         corr = 0,
                         perc_missing = 0,
                         X = NULL)
{
  # Whether or not a real design matrix is being used.
  real_design_matrix = !is.null(X)

  # Number of predictors:
  if(real_design_matrix)
  {
    p = ncol(X)
    n = nrow(X)
    stopifnot(p_0 < p)
  } else {
    p = 1 + p_1 + p_0
  }

  # Observation points:
  tau = seq(0, 1, length.out = m)

  # FLCs: orthonormalized polynomials
  F_true = cbind(1/sqrt(m),
                 poly(tau, K_true - 1))

  # Simulate the predictors:
  if(!real_design_matrix)
  {
    Xiid = matrix(rnorm(n = n*(p-1)), nr = n, nc = p-1)
    if(corr == 0){
      X = cbind(1,Xiid)
    } else {
      # Correlated predictors:
      ind_mat = matrix(rep(1:(p-1), p-1),nrow=p-1, byrow=FALSE);
      ind_diffs = abs(ind_mat - t(ind_mat))
      cov_mat = corr^ind_diffs
      # Cholesky:
      ch_cov_mat = chol(cov_mat)

      # Correlated predictors:
      X = cbind(1,
                t(crossprod(ch_cov_mat, t(Xiid))))
    }
  }

  # True coefficients:
  alpha_arr_true = array(0, c(K_true, p))

  # p = 1 Intercept coefficients: decaying importance
  alpha_arr_true[,1] = 1/(1:K_true)

  # Simulate the nonzero factors
  # Nonzero indices: if correlated predictors, or a real design matrix, space out the true ones
  nonzero_ind = 1:p_1; if(corr != 0 || real_design_matrix) nonzero_ind = round(seq(1, p-1, length.out = p_1))
  if(p_1 > 0){for(j in nonzero_ind){
    # Which factors are nonzero for predictor j?
    if(sparse_factors){ # Truncated Poisson(1)
      k_p = sample(1:K_true, truncdist::rtrunc(n = 1, spec = 'pois', a = 1, b = K_true, lambda = 1))
    } else k_p = 1:K_true

    # Simulate true values of the (nonzero) regression coefficients (decaying importance)
    alpha_arr_true[k_p, j+1] = rnorm(n = length(k_p), mean = 0, sd = 1/k_p)
  }}

  # True coefficient functions:
  alpha_tilde_true = F_true %*% alpha_arr_true # m x p

  # True factors: add Gaussian (subject-specific) noise (decaying importance)
  Beta_true = matrix(0, nrow = n, ncol = K_true)
  for(k in 1:K_true) Beta_true[,k] = X%*%alpha_arr_true[k,] + rnorm(n = n, sd = 1/k)

  # True FTS:
  Y_true = tcrossprod(Beta_true, F_true)

  # Noise SD based on RSNR:
  sigma_true = sd(Y_true)/RSNR

  # Observed data:
  Y = Y_true + sigma_true*rnorm(m*n)

  # Remove any observation points:
  if(perc_missing > 0 ) Y[sample(1:length(Y), perc_missing*length(Y))] = NA

  list(Y = Y, X = X, tau = tau,
       Y_true = Y_true, alpha_tilde_true = alpha_tilde_true,
       alpha_arr_true = alpha_arr_true, Beta_true = Beta_true, F_true = F_true, sigma_true = sigma_true)
}
#----------------------------------------------------------------------------
# Initialize the factors and FLCs using a SVD
  # Inputs: Y, tau, K (see fdlm() for details)
# Returns a list of the main parameters in the model:
# Beta, d, splineInfo (spline basis matrices), and the imputed data matrix Y0
fdlm_init_d = function(Y, tau, K){

  # Convert to matrix, if necessary:
  tau = as.matrix(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  # And the dimensions:
  n = nrow(Y); m = ncol(Y); d = ncol(tau)

  # Compute basic quantities for the FLC splines:
  splineInfo = getSplineInfo_d(tau = tau01,
                               m_eff = floor(median(rowSums(!is.na(Y)))),
                               orthonormalize = TRUE)

  # For initialization: impute (this is a bit crude)
  Y0 = matrix(NA, nr = n, nc = m) # Storage for imputed initialization data matrix
  allMissing.t = (rowSums(!is.na(Y))==0)   # Boolean indicator of times at which no points are observed
  # First: for all times at which we observe a curve, impute the full curve (across tau)
  #Y0[!allMissing.t,] = t(apply(Y[!allMissing.t,], 1, function(x) splinefun(tau01, x, method='natural')(tau01)))
  Y0[!allMissing.t,] = t(apply(Y[!allMissing.t,], 1, function(x) splinefun(seq(0, 1, length.out = m), x, method='natural')(seq(0, 1, length.out = m))))
  # Second: impute any times for which no curve is observed (i.e., impute across time)
  Y0 = apply(Y0, 2, function(x){splinefun(1:n, x, method='natural')(1:n)})

  # Compute SVD of the (completed) data matrix:
  # (Complete) data matrix, projected onto basis:
  YB0 = Y0%*%splineInfo$Bmat%*%chol2inv(chol(splineInfo$BtB))
  singVal = svd(YB0)

  # If K is unspecified, select based on cpv
  if(is.null(K)){
    # Cumulative sums of the s^2 proportions (More reasonable when Y has been centered)
    cpv = cumsum(singVal$d^2/sum(singVal$d^2))
    K = max(2, which(cpv >= 0.99)[1])
  }

  # Check to make sure K is less than the number of basis coefficients!
  if(K >= ncol(YB0)){
    warning(paste("K must be less than the number of basis functions used; reducing from K =",K, "to K =", ncol(YB0) - 1))
    K = ncol(YB0) - 1
  }

  # Basis coefficients of FLCs:
  Psi0 = as.matrix(singVal$v[,1:K])

  # Initial FLCs:
  F0 = splineInfo$Bmat%*%Psi0

  # Factors:
  Beta0 = as.matrix((singVal$u%*%diag(singVal$d))[,1:K])

  # Initialize all curves to have positive sums (especially nice for the intercept)
  negSumk = which(colSums(F0) < 0); Psi0[,negSumk] = -Psi0[,negSumk]; Beta0[,negSumk] = -Beta0[,negSumk]

  list(Beta = Beta0, Psi = Psi0, splineInfo = splineInfo, Y0 = Y0)
}
#----------------------------------------------------------------------------
# Construct the spline basis and penalty matrices
#
# Given input points in \code{d} dimensions, construct a low-rank thin plate spline basis matrix
# and penalty matrix.
#
# @param tau \code{m x d} matrix of coordinates, where \code{m} is the number of observation points and \code{d} is the dimension
# @param m_eff the effective number of observation points;
# (e.g., the median number of observation points when missing observations)
# @param orthonormalize logical; if TRUE, orthornomalize the basis matrix
#
# @note The knot locations are selected using a space-filling algorithm.
#
# @import fields
getSplineInfo_d = function(tau, m_eff = NULL, orthonormalize = TRUE){

  # Just in case, reform as matrix
  tau = as.matrix(tau)

  # Number of observation points
  m = nrow(tau)

  # Dimension:
  d = ncol(tau)

  # Order of derivative in penalty:
  m_deriv = 2

  # This is the linear component
  X = cbind(1, tau)

  # Number of effective observation points:
  if(is.null(m_eff)) m_eff = m

  # Number of knots: if more than 25 effective observation points, likely can use fewer knots
  if(m_eff > 25){
    # Guaranteed to be between 20 and 150 knots (but adjust as needed)
    num_knots = max(20, min(ceiling(m_eff/4), 150))
  } else num_knots = max(3, m_eff)

  # Locations of knots:
  if(num_knots < m){
    # Usual case: fewer knots than TOTAL observation points
    if(d == 1){
      # d = 1, just use quantiles of the observed data points:
      knots = as.matrix(quantile(unique(tau), seq(0,1,length=(num_knots+2))[-c(1,(num_knots+2))]))
    } else {
      # d > 1, use space-filling algorithm:
      knots = cover.design(tau, num_knots)$design
    }
  } else knots = tau

  # For the penalty matrix, need to compute distances between obs. points and knots
  dist.mat <- matrix(0, num_knots, num_knots); dist.mat[lower.tri(dist.mat)] <- dist(knots); dist.mat <- dist.mat + t(dist.mat)
  if(d%%2 == 0){
    # Even dim:
    Omega = dist.mat^(2*m_deriv - d)*log(dist.mat)
  } else {
    # Odd dim:
    Omega = dist.mat^(2*m_deriv - d)
  }
  # For numerical stability:
  diag(Omega) = 0

  # Compute the "random effects" matrix
  Zk = matrix(0, nrow=m, ncol=num_knots)
  for (k in 1:num_knots){
    di = sqrt(rowSums((tau - matrix(rep(knots[k,], each = m), nrow=m))^2)) # di = 0; for(j in 1:d) di = di + (tau[,j] - knots[k,j])^2; di = sqrt(di)
    if(d%%2 == 0){# Even dim:
      Zk[,k] = di^(2*m_deriv - d)*log(di)
    } else { # Odd dim:
      Zk[,k] = di^(2*m_deriv - d)
    }
  }
  Zk[is.nan(Zk)] = 0

  # Natural constraints, if necessary:
  if(num_knots > m - 1){Q2 = qr.Q(qr(X), complete=TRUE)[,-(1:2)]; Zk = Zk%*%Q2; Omega = crossprod(Q2, Omega)%*%Q2}

  # SVD of penalty matrix
  # So that the "random effects" have diagonal prior variance
  svd.Omega = svd(Omega)
  sqrt.Omega = t(svd.Omega$v %*%(t(svd.Omega$u)*sqrt(svd.Omega$d)))
  Z = t(solve(sqrt.Omega,t(Zk)))

  # Now combine the linear and nonlinear pieces to obtain the matrix of basis functions evaluated at the obs. points
  Bmat = cbind(X, Z);

  # The penalty matrix:
  Omega = diag(c(rep(0, ncol(X)), rep(1, ncol(Z))))

  if(orthonormalize){
    # QR decomposition:
    qrd = qr(Bmat, complete = TRUE);  R.t = t(qr.R(qrd));
    # Update hte basis and the penalty matrix:
    Bmat = qr.Q(qrd); Omega = forwardsolve(R.t, t(forwardsolve(R.t, Omega, upper.tri = FALSE)), upper.tri = FALSE)

    BtB = diag(1, ncol(Bmat))
  } else BtB = crossprod(Bmat)

  # Return the matrix, the penalty, and the cross product (of the basis)
  list(Bmat = Bmat, Omega = Omega, BtB = BtB)
}
#----------------------------------------------------------------------------
#' Compute the posterior distrubution for the regression coefficient functions
#'
#' Given posterior samples for the loading curves \code{fk} and the
#' regression coefficient factors \code{alpha_j} for a predictor \code{j},
#' compute the posterior distribution of the corresponding regression coefficient function.
#'
#' @param post_fk \code{Nsims x m x K} matrix of posterior draws of the loading curve matrix
#' @param post_alpha_j \code{Nsims x K} matrix of posterior draws of the regression coefficient factors
#'
#' @return \code{Nsims x m} matrix of posterior draws of the regression coefficient function
#' @export
get_post_alpha_tilde = function(post_fk, post_alpha_j){

  # Compute dimensions:
  Nsims = dim(post_fk)[1]; m = dim(post_fk)[2]

  # Store the result:
  post_alpha_tilde_j = matrix(0, nrow = Nsims, ncol = m)
  for(ni in 1:Nsims) post_alpha_tilde_j[ni,] = post_fk[ni,,]%*%post_alpha_j[ni,]

  post_alpha_tilde_j
}
#----------------------------------------------------------------------------
#' Compute Global Bayesian P-Values
#'
#' Given posterior samples for the loading curves \code{fk} and the
#' regression coefficient factors \code{alpha},
#' compute Global Bayesian P-Values for all regression coefficient functions
#'
#' @param post_fk \code{Nsims x m x K} matrix of posterior draws of the loading curve matrix
#' @param post_alpha_j \code{Nsims x p x K} matrix of posterior draws of the regression coefficient factors
#'
#' @return \code{p x 1} vector of Global Bayesian P-Values
#' @export
fosr_gbpv = function(post_fk, post_alpha){
  p = dim(post_alpha)[2]
  gbpv = numeric(p)
  for(j in 1:p){
    post_alpha_tilde_j = get_post_alpha_tilde(post_fk, post_alpha[,j,])
    gbpv[j] = min(simBaS(post_alpha_tilde_j))
  }
  gbpv
}
#----------------------------------------------------------------------------
# Compute Simultaneous Credible Bands
#
# Compute (1-alpha)\% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
#
# @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
# @param alpha confidence level
#
# @return \code{m x 2} matrix of credible bands; the first column is the lower band, the second is the upper band
#
# @note The input needs not be curves: the simultaneous credible "bands" may be computed
# for vectors. The resulting credible intervals will provide joint coverage at the (1-alpha)%
# level across all components of the vector.
#
credBands = function(sampFuns, alpha = .05){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # Compute the (1-alpha) sample quantile:
  Malpha = quantile(Maxfx, 1-alpha)

  # Finally, store the bands in a (m x 2) matrix of (lower, upper)
  cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx)
}
#----------------------------------------------------------------------------
# Compute Simultaneous Band Scores (SimBaS)
#
# Compute simultaneous band scores (SimBaS) from Meyer et al. (2015, Biometrics).
# SimBaS uses MC(MC) simulations of a function of interest to compute the minimum
# alpha such that the joint credible bands at the alpha level do not include zero.
# This quantity is computed for each grid point (or observation point) in the domain
# of the function.
#
# @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#
# @return \code{m x 1} vector of simBaS
#
# @note The input needs not be curves: the simBaS may be computed
# for vectors to achieve a multiplicity adjustment.
#
# @note The minimum of the returned value, \code{PsimBaS_t},
# over the domain \code{t} is the Global Bayesian P-Value (GBPV) for testing
# whether the function is zero everywhere.
#
simBaS = function(sampFuns){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # And now compute the SimBaS scores:
  PsimBaS_t = rowMeans(sapply(Maxfx, function(x) abs(Efx)/SDfx <= x))

  PsimBaS_t
}
#----------------------------------------------------------------------------
# Estimate the remaining time in the MCMC based on previous samples
# @param nsi Current iteration
# @param timer0 Initial timer value, returned from \code{proc.time()[3]}
# @param nsims Total number of simulations
# @param nrep Print the estimated time remaining every \code{nrep} iterations
# @return Table of summary statistics using the function \code{summary}
computeTimeRemaining = function(nsi, timer0, nsims, nrep=100){

  # Only print occasionally:
  if(nsi%%nrep == 0 || nsi==20) {
    # Current time:
    timer = proc.time()[3]

    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)

    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec

    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(round(secRemaining/60, 2), "minutes remaining"))
      } else print(paste(round(secRemaining, 2), "seconds remaining"))
    }
  }
}
#----------------------------------------------------------------------------
# Summarize of effective sample size
#
# Compute the summary statistics for the effective sample size (ESS) across
# posterior samples for possibly many variables
#
# @param postX An array of arbitrary dimension \code{(nsims x ... x ...)}, where \code{nsims} is the number of posterior samples
# @return Table of summary statistics using the function \code{summary()}.
#
# @examples
# # ESS for iid simulations:
# rand_iid = rnorm(n = 10^4)
# getEffSize(rand_iid)
#
# # ESS for several AR(1) simulations with coefficients 0.1, 0.2,...,0.9:
# rand_ar1 = sapply(seq(0.1, 0.9, by = 0.1), function(x) arima.sim(n = 10^4, list(ar = x)))
# getEffSize(rand_ar1)
#
# @import coda
getEffSize = function(postX) {
  if(is.null(dim(postX))) return(effectiveSize(postX))
  summary(effectiveSize(as.mcmc(array(postX, c(dim(postX)[1], prod(dim(postX)[-1]))))))
}
#----------------------------------------------------------------------------
# Compute the ergodic (running) mean.
# @param x vector for which to compute the running mean
# @return A vector \code{y} with each element defined by \code{y[i] = mean(x[1:i])}
# @examples
# # Compare:
# ergMean(1:10)
# mean(1:10)
#
## Running mean for iid N(5, 1) samples:
# x = rnorm(n = 10^4, mean = 5, sd = 1)
# plot(ergMean(x))
# abline(h=5)
ergMean = function(x) {cumsum(x)/(1:length(x))}
#----------------------------------------------------------------------------
# Compute the log-odds
# @param x scalar or vector in (0,1) for which to compute the (componentwise) log-odds
# @return A scalar or vector of log-odds
# @examples
# x = seq(0, 1, length.out = 10^3)
# plot(x, logit(x))
logit = function(x) {
  if(any(abs(x) > 1)) stop('x must be in (0,1)')
  log(x/(1-x))
}
#----------------------------------------------------------------------------
# Compute the inverse log-odds
# @param x scalar or vector for which to compute the (componentwise) inverse log-odds
# @return A scalar or vector of values in (0,1)
# @examples
# x = seq(-5, 5, length.out = 10^3)
# plot(x, invlogit(x))
invlogit = function(x) exp(x - log(1+exp(x))) # exp(x)/(1+exp(x))
#----------------------------------------------------------------------------
# Compute a block diagonal matrix w/ constant blocks
#
# The function returns kronecker(diag(nrep), Amat), but is computed more efficiently
# @param Amat matrix to populate the diagaonal blocks
# @param nrep number of blocks on the diagonal
#----------------------------------------------------------------------------
blockDiag = function(Amat, nrep){
  nr1 = nrow(Amat); nc1 = ncol(Amat)
  fullMat = matrix(0, nr = nr1*nrep, nc = nc1*nrep)
  rSeq = seq(1, nr1*nrep + nr1, by=nr1) # row sequence
  cSeq = seq(1, nc1*nrep + nc1, by=nc1) # col sequence
  for(i in 1:nrep) fullMat[rSeq[i]:(rSeq[i+1] - 1),  cSeq[i]:(cSeq[i+1] - 1)] = Amat

  fullMat
}
#----------------------------------------------------------------------------
#' Plot a curve given posterior samples
#'
#' Plot the posterior mean, simultaneous and pointwise 95\% credible bands
#' for a curve given draws from the posterior distribution
#' @param post_f \code{Ns x m} matrix of \code{Ns} posterior simulations
#' of the curve at \code{m} points
#' @param tau \code{m x 1} vector of observation points
#' @param alpha confidence level for the bands
#' @param include_joint logical; if TRUE, include joint bands (as well as pointwise)
#' @param main text for title plot
#' @export
plot_curve = function(post_f, tau = NULL, alpha = 0.05, include_joint = TRUE, main = "Posterior Mean and Credible Bands"){

  Ns = nrow(post_f); m = ncol(post_f)

  if(is.null(tau)) tau = 1:m

  par(mfrow = c(1, 1), mai = c(1, 1, 1, 1))

  # Pointwise intervals:
  dcip = dcib = HPDinterval(as.mcmc(post_f), prob = 1 - alpha);

  # Joint intervals, if necessary:
  if(include_joint) dcib = credBands(post_f, alpha = alpha)

  f_hat = colMeans(post_f)

  plot(tau, f_hat, type = "n", ylim = range(dcib, dcip, na.rm = TRUE),
       xlab = expression(tau), ylab = "", main = main,
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
  if(include_joint) polygon(c(tau, rev(tau)), c(dcib[, 2], rev(dcib[, 1])), col = "gray50",
                            border = NA)
  polygon(c(tau, rev(tau)), c(dcip[, 2], rev(dcip[, 1])), col = "grey",
          border = NA)
  lines(tau, f_hat, lwd = 8, col = "cyan")
}
#----------------------------------------------------------------------------
#' Plot the factors
#'
#' Plot posterior mean of the factors together with the simultaneous and pointwise
#' 95\% credible bands.
#'
#' @param post_beta the \code{Nsims x n x K} array of \code{Nsims} draws from the posterior
#' distribution of the \code{n x K} matrix of factors, \code{beta}
#' @param subj \code{n x 1} vector of subject IDs or labels
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics abline lines  par plot polygon
#' @import coda
#' @export
plot_factors = function(post_beta, subj = NULL){
  K = dim(post_beta)[3] # Number of factors
  if(is.null(subj)) subj = seq(0, 1, length.out = dim(post_beta)[2])

  dev.new(); par(mai = c(.8,.9,.4,.4), bg = 'gray90');
  plot(subj, post_beta[1,,1], ylim = range(post_beta), xlab = 'Subject', ylab = '', main = paste('Subject-Specific Factors', sep=''), type='n', cex.lab = 2, cex.axis=2,cex.main=2)
  abline(h = 0, lty=3, lwd=2);
  for(k in K:1){
    cb = credBands(as.mcmc(post_beta[,,k])); ci = HPDinterval(as.mcmc(post_beta[,,k]));
    polygon(c(subj, rev(subj)), c(cb[,2], rev(cb[,1])), col='grey50', border=NA);
    polygon(c(subj, rev(subj)), c(ci[,2], rev(ci[,1])), col='grey', border=NA);
    lines(subj,colMeans(post_beta[,,k]), lwd=8, col=k)
  }
}
#----------------------------------------------------------------------------
#' Plot the factor loading curves
#'
#' Plot posterior mean of the factor loading curves together with the simultaneous
#' and pointwise 95\% credible bands.
#'
#' @param post_fk the \code{Nsims x m x K} array of \code{Nsims} draws from the posterior
#' distribution of the \code{m x K} matrix of FLCs, \code{fk}
#' @param tau \code{m x 1} vector of observation points
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics abline lines  par plot polygon
#' @import coda
#' @export
plot_flc = function(post_fk, tau = NULL){
  K = dim(post_fk)[3] # Number of factors
  if(is.null(tau)) tau = seq(0, 1, length.out = dim(post_fk)[2])

  dev.new(); par(mai = c(.9,.9,.4,.4), bg = 'gray90');
  plot(tau, post_fk[1,,1], ylim = range(post_fk), xlab = expression(tau), ylab = '', main = 'Factor Loading Curves', type='n', cex.lab = 2, cex.axis=2,cex.main=2)
  abline(h = 0, lty=3, lwd=2);
  for(k in K:1){
    # Credible intervals:
    ci = HPDinterval(as.mcmc(post_fk[,,k]));
    # Credible bands (w/ error catch):
    cb = try(credBands(as.mcmc(post_fk[,,k])), silent = TRUE)
    if(class(cb) == "try-error") cb = ci
    polygon(c(tau, rev(tau)), c(cb[,2], rev(cb[,1])), col='grey50', border=NA);
    polygon(c(tau, rev(tau)), c(ci[,2], rev(ci[,1])), col='grey', border=NA);
    lines(tau,colMeans(post_fk[,,k]), lwd=8, col=k)
  }
}
#----------------------------------------------------------------------------
#' Plot the Bayesian curve fitted values
#'
#' Plot the curve posterior means with posterior credible intervals (pointwise and joint),
#' the observed data, and true curves (if known)
#'
#' @param y the \code{n x 1} vector of time series observations
#' @param mu the \code{n x 1} vector of fitted values, i.e., posterior expectation of the mean
#' @param postY the \code{nsims x n} matrix of posterior draws from which to compute intervals
#' @param y_true the \code{n x 1} vector of points along the true curve
#' @param t01 the observation points; if NULL, assume \code{n} equally spaced points from 0 to 1
#' @param include_joint_bands logical; if TRUE, compute simultaneous credible bands
#'
#' @examples
#' # FIXME
#'
#' @import coda
#' @export
plot_fitted = function(y, mu, postY, y_true = NULL, t01 = NULL, include_joint_bands = FALSE){

  # Time series:
  n = length(y);
  if(is.null(t01)) t01 = seq(0, 1, length.out=n)

  # Credible intervals/bands:
  #dcip = HPDinterval(as.mcmc(postY)); dcib = credBands(postY)
  dcip = dcib = t(apply(postY, 2, quantile, c(0.05/2, 1 - 0.05/2)));
  if(include_joint_bands) dcib = credBands(postY)

  # Plot
  dev.new(); par(mfrow=c(1,1), mai = c(1,1,1,1))
  plot(t01, y, type='n', ylim=range(dcib, y, na.rm=TRUE), xlab = 't', ylab=expression(paste("y"[t])), main = 'Fitted Values: Conditional Expectation', cex.lab = 2, cex.main = 2, cex.axis = 2)
  polygon(c(t01, rev(t01)), c(dcib[,2], rev(dcib[,1])), col='gray50', border=NA)
  polygon(c(t01, rev(t01)), c(dcip[,2], rev(dcip[,1])), col='grey', border=NA)
  if(!is.null(y_true))  lines(t01, y_true, lwd=8, col='black', lty=6);
  lines(t01, y, type='p');
  lines(t01, mu, lwd=8, col = 'cyan');
}
#----------------------------------------------------------------------------
# Univariate Slice Sampler from Neal (2008)
#
# Compute a draw from a univariate distribution using the code provided by
# Radford M. Neal. The documentation below is also reproduced from Neal (2008).
#
# @param x0    Initial point
# @param g     Function returning the log of the probability density (plus constant)
# @param w     Size of the steps for creating interval (default 1)
# @param m     Limit on steps (default infinite)
# @param lower Lower bound on support of the distribution (default -Inf)
# @param upper Upper bound on support of the distribution (default +Inf)
# @param gx0   Value of g(x0), if known (default is not known)
#
# @return  The point sampled, with its log density attached as an attribute.
#
# @note The log density function may return -Inf for points outside the support
# of the distribution.  If a lower and/or upper bound is specified for the
# support, the log density function will not be called outside such limits.
uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # Check the validity of the arguments.

  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }

  # Keep track of the number of calls made to this function.
  #uni.slice.calls <<- uni.slice.calls + 1

  # Find the log density at the initial point, if not already known.

  if (is.null(gx0))
  { #uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0)
  }

  # Determine the slice level, in log terms.

  logy <- gx0 - rexp(1)

  # Find the initial interval to sample from.

  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }

    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }

  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J

    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }

    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }

  # Shrink interval to lower and upper bounds.

  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }

  # Sample from the interval, shrinking it on each rejection.

  repeat
  {
    x1 <- runif(1,L,R)

    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)

    if (gx1>=logy) break

    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }

  # Return the point sampled, with its log density attached as an attribute.

  attr(x1,"log.density") <- gx1
  return (x1)

}

# Just add these for general use:
#' @importFrom stats quantile rgamma rnorm sd splinefun var rexp runif
NULL
