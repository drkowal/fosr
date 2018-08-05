#----------------------------------------------------------------------------
#' Decoupling shrinkage and selection for function-on-scalars regression
#'
#' For a functional response and scalar predictors, construct a posterior
#' summary that balances predictive accuracy and sparsity. Given posterior
#' draws of regression coefficients (or coefficient functions) from a FOSR model,
#' use a suitably-defined loss function to select important variables for prediction.
#'
#'
#' @param X \code{n x p} matrix of predictors
#' @param post_alpha \code{Nsims x p x K} array of \code{Nsims} posterior draws
#' of the \code{p} predictors for each of \code{K} factors
#' @param post_trace_sigma_2 \code{Nsims x 1} vector of posterior draws of the trace
#' of the (marginal) covariance (see below for details)
#' @param weighted logical; if TRUE, use weighted group lasso (recommended)
#' @param alpha_level coverage for the credible interval on the proportion of
#' variance explained
#' @param remove_int logical; if TRUE, remove the intercept term from model comparisons
#' @param include_plot logical; if TRUE, include a plot showing proportion of variability
#' explained against model size
#' @param include_model_list; if TRUE, include model_list in return--a boolean matrix
#' of models of different sizes suggested by DSS
#'
#' @note This function is value for the regression functions (m-dimensional) as well as the
#' regression factors (K-dimensional). Since K << m, the latter is much faster.
#'
#' @note The matrix of predictors, \code{X}, may be different from the given matrix
#' in the data; i.e., we may have a different set of design points for prediction.
#'
#' @note \code{post_trace_sigma_2} is the (posterior samples of)
#' the trace of the error covariance matrix jointly across subjects i=1,...,n
#' and observations j=1,...,m, after marginalizing out the random effects \code{gamma_ik}.
#' This is given by \code{nm x sigma_e^2} + \code{sum_ik sigma_gamma_ik^2},
#' where the second term is necessary only when random effects are included in the model
#' AND integrated over in the predictive distribution.
#'
#' @return alpha_dss a \code{p x K} matrix of (sparse) regression coefficents; if
#' include_model_list is TRUE, return a list of alpha_dss and model_list, a boolean matrix
#' of possible different sized models
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
#' # Run the DSS:
#' alpha_dss = fosr_select(X = X,
#'                        post_alpha = out$alpha,
#'                        post_trace_sigma_2 = n*m*out$sigma_e^2 + apply(out$sigma_g^2, 1, sum))
#' # Variables selected:
#' (select_dss = which(apply(alpha_dss, 1, function(x) any(x != 0))))
#'
#' @import gglasso lars
#' @export
fosr_select = function(X,
                       post_alpha,
                       post_trace_sigma_2,
                       weighted = TRUE,
                       alpha_level = 0.10,
                       remove_int = TRUE,
                       include_plot = TRUE,
                       include_model_list = FALSE){
  #----------------------------------------------------------------------------
  # Run some checks:
  if(ncol(X)!= dim(post_alpha)[2])
    stop('Dimensions of X and post_alpha are not compatible')

  if(length(post_trace_sigma_2) != dim(post_alpha)[1])
    stop('post_trace_sigma_2 and post_alpha must have the same number of simulations')

  if(remove_int && length(unique(X[,1])) > 1)
    stop('To remove an intercept (remove_int = TRUE), the first column of X must be constant (e.g., all ones)')
  #----------------------------------------------------------------------------

  # Posterior mean:
  alpha_hat = colMeans(post_alpha)

  #----------------------------------------------------------------------------
  # Group lasso variable selection:

  # Weights:
  w = NULL; if(weighted) w = 1/rowSums(alpha_hat^2)

  # (Adaptive) group lasso:
  alpha_lam = dss_select(Y = X%*%alpha_hat, X = X, w = w)

  #----------------------------------------------------------------------------
  # Compute the variability explained

  # Remove intercept for variability explained?
  if(remove_int) {
    alpha_int = alpha_lam[,1,] # Save the intercept terms (along the lambda-path)
    alpha_lam = alpha_lam[,-1,];
    X = X[,-1];
    post_alpha = post_alpha[,-1,]
  }

  # Compute the proportion of variability explained part:
  rhoList = computeRho(beta_path = alpha_lam,
                       XX = X,
                       post_beta = post_alpha,
                       post_trace_sigma_2 = post_trace_sigma_2)
  rho_lam0 = rhoList$rho_lam0; rho_lam2 = rhoList$rho_lam2;

  # Number of lambda values in the path:
  L = ncol(rho_lam2)

  #----------------------------------------------------------------------------
  # Select the variables

  # HPD interval for variability explained in the sparsified models:
  ci_2 = HPDinterval(as.mcmc(rho_lam2), prob = 1 - alpha_level)

  # Cutoff: use posterior mean for the full model
  rho_lam0_hat = mean(rho_lam0)

  # Identify the models in the lambda-path for which the HPD interval
    # contains the posterior mean of the full model:
  if(!any(ci_2[,2] >rho_lam0_hat)){
    warning('Credible intervals for rho_lam do not contain the posterior mean of rho_0;
            using the posterior mean of rho_lam corresponding to
            the largest model in the lambda-path instead.')
    mrange = which(ci_2[,2] > mean(rho_lam2[,L])) # Model range
  } else mrange = which(ci_2[,2] > rho_lam0_hat) # Model range

  # Count the number of variables (for each lambda) that are nonzero for ALL k=1,...,K:
  model_size = apply(alpha_lam, 3, function(a) sum(colSums(a) != 0))
  #model_size = numeric(L); for(ell in 1:L) model_size[ell] = sum(alpha_lam[1,,ell] != 0) #sum(colSums(alpha_lam[,,ell]) != 0)

  # Select the sparsest model that is "equivalent" to the full model
    # in this HPD interval overlap sense
  ind = max(which(model_size[mrange] == min(model_size[mrange])))
  alpha_dss = t(alpha_lam[,,mrange[ind]])

  # Include the intercept:
  if(remove_int) alpha_dss = rbind(alpha_int[,mrange[ind]],
                                   alpha_dss)

  # Number of variables selected:
  p_dss = sum(rowSums(alpha_dss) != 0)
  print(paste('DSS Selected', p_dss, 'variables (including the intercept)'))

  #----------------------------------------------------------------------------
  # Include model list:
  if(include_model_list){
    # pick the model sizes we care about
    idxs = which(c(model_size, max(model_size)) - c(0, model_size) != 0)

    # for each of the model sizes, record the selected model
    # The dimension is number of models we are considering x (p_0 + p_1)
    model_list = array(FALSE, dim = c(length(idxs), dim(alpha_lam)[2]))
    rownames(model_list) = idxs
    for(idx in idxs){
      model = which(colSums(alpha_lam[,,idx]) != 0)
      model_list[as.character(idx),model] = TRUE
    }
    rownames(model_list) = NULL

    # Make sure to add back the intercept if needed
    # remove_int means remove intercept for variance calculated--
    #  which means the intercept is to be included in the model
    if(remove_int){
      model_list = cbind(TRUE, model_list)
    }
  }

  #----------------------------------------------------------------------------
  # Include a summary plot:
  if(include_plot){
    # HPD interval for full model:
    ci_0 = HPDinterval(as.mcmc(rho_lam0), prob = 1 - alpha_level)

    plot(model_size, model_size, ylim = range(0, 1), type='n',
         xlab = 'Number of Predictors', ylab = expression(rho[lambda]^2))
    polygon(c(model_size, rev(model_size)), c(rep(ci_0[2], L), rev(rep(ci_0[1],L))), col = "grey",
            border = NA)
    abline(h = rho_lam0_hat, lty = 2, lwd=2)

    msize = unique(model_size)
    for(mi in 1:length(msize)){
      ind.mi = max(which((model_size == msize[mi]))) # match(msize[mi], model_size)
      lines(rep(msize[mi], 2), c(ci_2[ind.mi,1], ci_2[ind.mi,2]), lwd=2)
      #arrows(msize[mi], ci_2[ind.mi,1], msize[mi], ci_2[ind.mi,2], length=0.05, angle=90, code=3)
      lines(msize[mi], mean(rho_lam2[,ind.mi]), type = 'p', lwd=2)
    }
  }

  if(include_model_list){
    return(list(
      alpha_dss = alpha_dss,
      model_list = model_list))
  } else
    return(alpha_dss)
}
#----------------------------------------------------------------------------
# Group lasso for matrix regression
#
# Given a \code{N x M} data matrix \code{Y}
# and a \code{N x P} matrix of predictor \code{X},
# estimate the model
# \code{Y = XB + E}
# for \code{P x M} regression coefficient matrix \code{B}.
# The penalty applies a row-wise group lasso penalty, i.e.,
# for row \code{p}, the \code{M} elements \code{X[p,]} are regularized toward zero.
#
# @param Y \code{N x M} matrix of observations
# @param X \code{N x P} matrix of predictors
# @param w \code{P x 1} vector of weights; if NULL, use sqrt(M)
#
# @note The design matrix \code{X} may include an intercept, which will
# be left unpenalized.
#
# @return The solution path for \code{L} values of \code{lambda} in the form
# of an array of dimension \code{M x P x L}
# @import gglasso
dss_select = function(Y, X, w = NULL){

  # Number of observations:
  N = nrow(Y)

  # Number of replicates:
  M = ncol(Y)

  # Number of predictors (groups):
  P = ncol(X)

  # Construct the relevant terms for gglasso:
  yy = matrix(t(Y)) # vectorized
  XX = kronecker(X, diag(M))
  gx = rep(1:P, each = M)

  # Default choice for weights:
  if(is.null(w)) w = rep(sqrt(M), P)

  # Check: if there's an intercept, set the weights to zero
  w[which(apply(X, 2, function(x) length(unique(x)) == 1))] = 0

  # OLD CODE: Rescale XX if we use a weighted group lasso:
  #W = matrix(rep(1/w, times = M), ncol = M)
  #if(!is.null(W)) XX = t(as.numeric(t(W))*t(XX)) # XX%*%diag(as.numeric(t(W)))

  # Run the (adaptive) group lasso:
  beta_lam = as.matrix(gglasso(x = XX,
                               y = yy,
                               group = gx,
                               pf = w,
                               loss = "ls",
                               lambda.factor = 0, # 0.0001,
                               intercept = FALSE)$beta) # intercept should be in X
  L = ncol(beta_lam)

  # M x P x L (=#lambdas)
  beta_lam = array(beta_lam, c(M, P, L))

  # OLD CODE: Rescale beta if we use a weighted group lasso:
  #if(!is.null(W)){for(ell in 1:L) beta_lam[,,ell] = beta_lam[,,ell]*t(W)}

  return(beta_lam)
}
#----------------------------------------------------------------------------
# Lasso for matrix regression
#
# Given a \code{N x M} data matrix \code{Y}
# and a \code{N x P} matrix of predictor \code{X},
# estimate the model
# \code{Y = XB + E}
# for \code{P x M} regression coefficient matrix \code{B}.
# The penalty applies a lasso penalty to all elements of \code{B}.
#
# @param Y \code{N x M} matrix of observations
# @param X \code{N x P} matrix of predictors
# @param W \code{P x M} matrix of weights
#
# @return The solution path for \code{L} values of \code{lambda} in the form
# of an array of dimension \code{M x P x L}
# @import lars
dss_select_lasso = function(Y, X, W = NULL){

  # Number of observations:
  N = nrow(Y)

  # Number of replicates:
  M = ncol(Y)

  # Number of predictors (groups):
  P = ncol(X)

  # Construct the relevant terms for gglasso:
  yy = matrix(Y) # vectorized
  XX = blockDiag(X, M) # kronecker(diag(M), X)

  # Rescale XX if we use a weighted group lasso:
  if(!is.null(W)) XX = t(as.numeric(W)*t(XX)) # XX%*%diag(as.numeric(W))

  # Run the lasso:
  #beta_lam = as.matrix(glmnet(x = XX, y = yy, alpha = 1, intercept = TRUE, standardize = FALSE)$beta)
  beta_lam = t(as.matrix(lars(x = XX,
                              y = yy,
                              type = 'lasso',
                              intercept = TRUE,
                              normalize = FALSE,
                              use.Gram = I(P*M < 500))$beta))
  L = ncol(beta_lam)

  # M x P x L (=#lambdas)
  beta_lam = aperm(array(beta_lam, c(P, M, L)), c(2,1,3))

  # Rescale beta if we use a weighted group lasso:
  if(!is.null(W)){for(ell in 1:L) beta_lam[,,ell] = beta_lam[,,ell]*t(W)}

  return(beta_lam)
}
# Compute posterior predictive distribution of \code{rho.lambda}
#
# Using the (group) lasso solution path, compute the proportion of
# variability explained by the sparsified models (relative to the full model)
# for each MCMC simulation.
#
# @param beta_path \code{M x P x L} array of regression coefficients \code{beta}
# along the solution path of length \code{L}
# @param XX \code{N x P} matrix of predictors
# @param post_beta \code{Nsims x P x M} array of posterior draws of \code{beta}
# @param post_trace_sigma_2 \code{Nsims x 1} vector of posterior draws of the trace
# of the (marginal) covariance (see below for details)
#
# @note \code{post_trace_sigma_2} is the (posterior samples of)
# the trace of the error covariance matrix jointly across subjects i=1,...,n
# and observations j=1,...,m, after marginalizing out the random effects \code{gamma_ik}.
# This is given by \code{nm x sigma_e^2} + \code{sum_ik sigma_gamma_ik^2},
# where the second term is necessary only when random effects are included in the model
# AND integrated over in the predictive distribution.
#
# @return A list containing \code{rho_lam0} and \code{rho_lam2},
# corrosponding to rho2 for the full model and the sparsified model
# (for each value of \code{lambda} in the solution path).
computeRho = function(beta_path, XX, post_beta, post_trace_sigma_2){
  # Use the same notation as before:
  N = nrow(XX)            # Number of observations
  M = dim(beta_path)[1];  # Number of factors (rows)
  P = dim(beta_path)[2];  # Number of predictors (excluding the intercept)
  L = dim(beta_path)[3];  # Number of lambda values considered
  Ns = length(post_trace_sigma_2) # Number of MCMC simulations

  # Posterior mean:
  beta_hat = colMeans(post_beta)

  # Construct storage matrices:
  rho_lam2 = matrix(0, nrow = Ns, ncol = L)
  rho_lam0 = numeric(Ns)

  for(ni in 1:Ns){
    XB = XX%*%post_beta[ni,,]

    # Squared norm:
    XB2 = sum((XB)^2)

    # Zeroth version:
    rho_lam0[ni] = XB2/(XB2 + post_trace_sigma_2[ni] + sum((XB - XX%*%beta_hat)^2))

    # ellth version:
    #for(ell in 1:L) rho_lam2[ni, ell] = XB2/(XB2 + post_trace_sigma_2[ni] + sum((XB - tcrossprod(XX, beta_path[,,ell]))^2))
    rho_lam2[ni, ] = XB2/(XB2 + post_trace_sigma_2[ni] + apply(beta_path, 3, function(b) sum((XB - tcrossprod(XX, b))^2)))
  }
  return(list(rho_lam0 = rho_lam0,
              rho_lam2 = rho_lam2))
}
