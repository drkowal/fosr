#----------------------------------------------------------------------------
# Sample a Gaussian vector using the fast sampler of BHATTACHARYA et al.
#
# Sample from N(mu, Sigma) where Sigma = solve(crossprod(Phi) + solve(D))
# and mu = Sigma*crossprod(Phi, alpha):
#
# @param Phi \code{n x p} matrix (of predictors)
# @param Ddiag \code{p x 1} vector of diagonal components (of prior variance)
# @param alpha \code{n x 1} vector (of data, scaled by variance)
# @return Draw from N(mu, Sigma), which is \code{p x 1}, and is computed in \code{O(n^2*p)}
# @note Assumes D is diagonal, but extensions are available
# @keywords internal
sampleFastGaussian = function(Phi, Ddiag, alpha){

  # Dimensions:
  Phi = as.matrix(Phi); n = nrow(Phi); p = ncol(Phi)

  # Step 1:
  u = rnorm(n = p, mean = 0, sd = sqrt(Ddiag))
  delta = rnorm(n = n, mean = 0, sd = 1)

  # Step 2:
  v = Phi%*%u + delta

  # Step 3:
  w = solve(crossprod(sqrt(Ddiag)*t(Phi)) + diag(n), #Phi%*%diag(Ddiag)%*%t(Phi) + diag(n)
            alpha - v)

  # Step 4:
  theta =  u + Ddiag*crossprod(Phi, w)

  # Return theta:
  theta
}
#----------------------------------------------------------------------------
# Factor Loading Curve Sampling Algorithm
#
# Sample the factor loading curve basis coefficients subject to an orthonormality constraint.
# Additional linear constraints may be included as well.
#
# @param BtY \code{Lm x n} matrix \code{B.t()*Y} for basis matrix B
# @param Beta \code{n x K} matrix of factors
# @param Psi \code{Lm x K} matrix of previous factor loading curve coefficients
# @param BtB \code{Lm x Lm} matrix of \code{B.t()*B}
# @param Omega \code{Lm x Lm} prior precision/penalty matrix
# @param BtCon (optional) \code{Lm x Jc} matrix of additional constraints, pre-multiplied by B.t()
# @param lambda \code{K}-dimensional vector of prior precisions
# @param sigmat2 \code{n}-dimensional vector of time-dependent observation error variances
# @return Psi \code{Lm x K} matrix of factor loading curve coefficients
#
# @note This is a wrapper for Rcpp functions for the special cases of
# \code{K = 1} and whether or not additional (linear) constraints are included,
# i.e., whether or not \code{BtCon} is non-\code{NULL}.
# @keywords internal
fdlm_flc = function(BtY, Beta, Psi, BtB, Omega, BtCon = NULL, lambda, sigmat2){

  # Obtain the dimensions, in order of appearance:
  Lm = nrow(BtY); n = ncol(BtY); K = ncol(Beta);

  # Allow for scalar variance input
  if(length(sigmat2) == 1) sigmat2 = rep(sigmat2, n)

  # Check dimensions:
  if( (nrow(Beta) != n) ||
      (nrow(Psi) != Lm) || (ncol(Psi) != K) ||
      (nrow(BtB) != Lm) || (ncol(BtB) != Lm) ||
      (nrow(Omega) != Lm) || (ncol(Omega) != Lm) ||
      (length(lambda) != K) ||
      (length(sigmat2) != n)
  ) stop("Mismatched dimensions in FLC sampler")

  # No additional constraints (besides orthonormality of FLCs themselves)
  if(is.null(BtCon)){

    if(K == 1){
      # Special case: (FLC) orthogonality not necessary
      Psi = sampleFLC_1(BtY = BtY, Beta = Beta, Psi = Psi, BtB = BtB, Omega = Omega, lambda = lambda, sigmat2 = sigmat2)
    } else {
      Psi = sampleFLC(BtY = BtY, Beta = Beta, Psi = Psi, BtB = BtB, Omega = Omega, lambda = lambda, sigmat2 = sigmat2)
    }

  } else {
    # Additional constraints: orthogonal to BtCon


    # Special case: (FLC) orthogonality not necessary
    if(K == 1){
      Psi = sampleFLC_cons_1(BtY = BtY, Beta = Beta, Psi = Psi, BtB = BtB, Omega = Omega, BtCon = BtCon, lambda = lambda, sigmat2 = sigmat2)
    } else {
      Psi = sampleFLC_cons(BtY = BtY, Beta = Beta, Psi = Psi, BtB = BtB, Omega = Omega, BtCon = BtCon, lambda = lambda, sigmat2 = sigmat2)
    }

  }
}
#----------------------------------------------------------------------------
# Factor loading curve smoothing parameter sampler
#
# Sample the smoothing parameters for each factor loading curve.
#
# @param lambda \code{K}-dimensional vector of smoothing parameters (prior precisions)
# from previous MCMC iteration
# @param Psi \code{Lm x K} matrix of basis coefficients, where \code{Lm} is the number of
# basis functions and \code{K} is the number of factors
# @param Omega \code{Lm x Lm} penalty matrix; if NULL, assume it is diag(0, 0, 1,...,1)
# @param d dimension of \code{tau}; default is 1
# @param uniformPrior logical; when TRUE, use a uniform prior on prior standard deviations,
# \code{1/sqrt{lambda[k]}}; otherwise use independent Gamma(0.001, 0.001) prior for each \code{lambda[k]}
# @param orderLambdas logical; when TRUE, enforce the ordering constraint \code{lambda[1] > ... > lambda[K]}
# for identifiability

# @return The \code{K}-dimensional vector of samoothing parameters, \code{lambda}.
# @import truncdist
sample_lambda = function(lambda, Psi, Omega = NULL, d = 1, uniformPrior = TRUE, orderLambdas = TRUE){
  Lm = nrow(Psi); K = ncol(Psi)

  if(uniformPrior){shape0 = (Lm - d + 1 + 1)/2} else shape0 = (Lm - d - 1)/2 + 0.001; # for Gamma(0.001, 0.001) prior
  #if(uniformPrior){shape0 = (Lm + 1)/2} else shape0 = (Lm - 2)/2 + 0.001; # for Gamma(0.001, 0.001) prior

  for(k in 1:K){
    if(is.null(Omega)){rate0 = crossprod(Psi[-(1:(d+1)),k])/2} else rate0 = crossprod(Psi[,k], Omega)%*%Psi[,k]/2
    #if(is.null(Omega)){rate0 = crossprod(Psi[-(1:2),k])/2} else rate0 = crossprod(Psi[,k], Omega)%*%Psi[,k]/2

    if(!uniformPrior) rate0 = rate0 + 0.001  # for Gamma(0.001, 0.001) prior

    # Lower and upper bounds, w/ ordering constraints (if specified):
    if(orderLambdas){
      lam.l = 10^-8; lam.u = Inf; if(k != 1) lam.u = lambda[k-1];  # if(k != K) lam.l = lambda[k+1];
      lambda[k] = truncdist::rtrunc(1, 'gamma', a=lam.l, b=lam.u, shape=shape0, rate=rate0) # more stable, possibly faster
    } else lambda[k] = rgamma(1, shape = shape0, rate = rate0)
  }
  lambda
}
#--------------------------------------------------------------
# Multiplicative Gamma Process (MGP) Sampler
#
# Sample the global parameters, delta.h, with tau.h = cumprod(delta.h), from
# the MGP prior for factor models
#
#
# @param theta.jh the \code{p x K} matrix with entries theta.jh ~ N(0, tau.h^-1), j=1:p, h=1:K
# @param delta.h the \code{K}-dimensional vector of previous delta.h values, tau.h[h] = prod(delta.h[1:h])
# @param a1 the prior parameter for factor 1: delta.h[1] ~ Gamma(a1, 1)
# @param a2 the prior parameter for factors 2:K: delta.h[h] ~ Gamma(a2, 1) for h = 2:K
# @return \code{delta.h}, the \code{K}-dimensional vector of multplicative components,
# where tau.h[h] = prod(delta.h[1:h])
#
# @note The default \code{a1 = 2} and \code{a2 = 3} appears to offer the best performance
# in Durante (2017).
sampleMGP = function(theta.jh, delta.h, a1 = 2, a2 = 3){

  # Just in case:
  theta.jh = as.matrix(theta.jh)

  # Store the dimensions locally
  p = nrow(theta.jh); K = ncol(theta.jh)

  # Sum over the (squared) replicates:
  sum.theta.l = colSums(theta.jh^2)

  # h = 1 case is separate:
  tau.not.1 = cumprod(delta.h)/delta.h[1]
  delta.h[1] = rgamma(n = 1, shape = a1 + p*K/2,
                      rate = 1 + 1/2*sum(tau.not.1*sum.theta.l))
  # h > 1:
  if(K > 1){for(h in 2:K){
    tau.not.h = cumprod(delta.h)/delta.h[h]
    delta.h[h] = rgamma(n = 1, shape = a2 + p/2*(K - h + 1),
                        rate = 1 + 1/2*sum(tau.not.h[h:K]*sum.theta.l[h:K]))
  }}
  delta.h
}
