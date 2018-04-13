get_FLC = function(tau, nfunctions)
{
  nc = nfunctions
  nr = length(tau)

  # FLCs, orthogonalized, are polynomials centered at 1/2
  F_true = matrix(0, nrow = nr, ncol = nc)
  for(k in 1:nc)
    F_true[,k] = (tau - 0.5)^(k-1)

  # Orthonormalize:
  F_true = qr.Q(qr(F_true))

  # Positive:
  for(k in 1:nc)
  {
    if(sum(F_true[,k]) < 0)
      F_true[,k] = -F_true[,k]
  }

  return(F_true)
}

#----------------------------------------------------------------------------
# Simulate some data:
#----------------------------------------------------------------------------
gen_data_orig = function(
  nT = 200,    # number of time points (or subjsects)
  nM = 100,    # Number of observation points
  RSNR = 5,    # Root signal-to-noise ratio
  ar1 = 0,     # Autoregressive coefficient of the predictors X_j (i.e., autocorrelated?)
  K_true = 4,  # Number of factors
  p_0 = 2,     # Number of zero coefficients
  p_1 = 5)     # Number of nonzero "coefficients"
               #   (predictors j for which \alpha_{k,j} != 0 for some k)
{
  # Number of predictors:
  p = 1 + p_1 + p_0

  # Observation points:
  tau = seq(0, 1, length.out = nM)

  F_true = get_FLC(tau, nK)
  #plot(as.ts(F_true))

  # Simulate the predictors:
  if(ar1 == 0)
  {
    X = cbind(1,matrix(rnorm(n = nT*(p-1)), nr = nT, nc = p-1))
  }
  else
  {
    X = cbind(
      1,
      apply(
        matrix(0, nr = nT, nc = p-1),
        2,
        function(x)
        {
          arima.sim(n = nT, list(ar = ar1), sd = sqrt(1-ar1^2))
        }))
  }
  #plot(as.ts(X))

  # True coefficients:
  alpha.arr.true = array(0, c(p, K_true))

  # p = 1 Intercept coefficients: decaying importance
  alpha.arr.true[1,] = 1/(1:K_true)

  # Randomly select which k's are nonzero (k_p) for each predictor:
  if(p_1 > 0)
  {
    for(j in 1:p_1)
    {
      # Which factors are nonzero for predictor j? Truncated Poisson(1)
      k_p = sample(
        1:K_true,
        rtrunc(n = 1, spec = 'pois', b = K_true, lambda = 1))

      # Simulate true values of the (nonzero)
      #  regression coefficients (decaying importance)
      alpha.arr.true[j+1, k_p] =
        rnorm(n = length(k_p), mean = 0, sd = 1/k_p)
    }
  }

  # True factors: add Gaussian (subject-specific) noise (decaying importance)
  Beta_true = matrix(0, nrow = nT, ncol = K_true)
  for(k in 1:K_true)
  {
    Beta_true[,k] = X%*%alpha.arr.true[,k] + rnorm(n = nT, sd = 1/k)
  }

  # True FTS:
  Y_true = tcrossprod(Beta_true, F_true)

  # Noise SD based on RSNR:
  sigma_true = sd(Y_true)/RSNR

  # Observed data:
  Y = Y_true + sigma_true*rnorm(nM*nT)

  return(list(
    Y = Y,
    tau = tau,
    X = X,
    Y_true = Y_true,
    sigma_true = sigma_true,
    F_true = F_true,
    Beta_true = Beta_true,
    K_true = K_true,
    alpha_true = alpha.arr.true))
}

gen_data_group_predictors = function(
  nT = 200,    # number of time points (or subjsects)
  nM = 100,    # Number of observation points
  RSNR = 5,    # Root signal-to-noise ratio
  group_sizes = c(2, 2, 2), # the number of predictors in each group
  n_true_predictor_groups = 1,   # the number of predictor groups that are predictive
  sig_within_group = 0.3,        # variance of the groups
  n_loading_curves = length(group_sizes)) # number of loading curves to
                                          # generate the data from
{
  # Observation points:
  tau = seq(0, 1, length.out = nM)

  # add the intercept
  group_sizes = c(1, group_sizes)
  nG = length(group_sizes)
  nP = sum(group_sizes)
  nK = n_loading_curves
  nTG = n_true_predictor_groups + 1

  # get the indices to access the groups
  right_idx = cumsum(group_sizes)
  left_idx = c(1, cumsum(group_sizes) + 1)[1:nG]
  get_group_idx = function(g){ left_idx[g]:right_idx[g] }

  # simulate the predictors
  X = cbind(1, array(dim = c(nT, nP-1)))

  X_true = cbind(1, array(rnorm(nT*(nG-1)), c(nT, nG-1)))

  if(length(sig_within_group) == 1)
  {
    sig_within_group = c(0, rep(sig_within_group, nG-1))
  }
  else
  {
    stopifnot(length(sig_within_group) == nG-1)
  }

  for(g_idx in 1:nG)
  {
    X[,get_group_idx(g_idx)] =
      X_true[,g_idx] +
      rnorm(nT*group_sizes[g_idx], sd = sig_within_group[g_idx])
  }

  B = get_FLC(tau, nK)
  A = array(rnorm(nK*nTG), dim = c(nTG, nK))

  Y_true = X_true[,1:nTG] %*% tcrossprod(A, B)

  # Noise SD based on RSNR:
  sigma_true = sd(Y_true)/RSNR

  # Observed data:
  Y = Y_true + sigma_true*array(rnorm(nM*nT), c(nT, nM))

  return(list(
    Y = Y,
    tau = tau,
    X = X,
    X_true = X_true,
    Y_true = Y_true,
    sigma_true = sigma_true))
}

# how did they choose which functions are zero?!?!?!
gen_data_chen_goldsmith_16 = function(
  nT = 100,
  nP = 20, # nP_true has to be 3 since those beta(t) are fixed
  nD = 25,
  cov_mat = NULL) # they use G + I
{
  # rename variables saved in sysdata.rda
  beta.true = chen_beta.true
  cov.ranef = chen_cov.ranef
  cov.resid = chen_cov.resid

  if(is.null(cov_mat))
  {
    cov_mat = cov.resid + diag(nD)
  }

  beta.true = beta.true / 10

  xs = seq(0, 1, length.out = 25)
  beta1 = splinefun(xs, beta.true[1,])
  beta2 = splinefun(xs, beta.true[2,])
  beta3 = splinefun(xs, beta.true[3,])

  #beta1 = splinefun(
  #  seq(0, 1, 0.125),
  #  c(0.0125, 0.02, 0.022, 0.025, 0.025, 0, -.02, -.03, -.02))
  #beta2 = splinefun(
  #  seq(0, 1, 0.125),
  #  c(0.025, 0.03, 0.035, 0.035, 0.025, -0.025, -0.09, -.12, -.115))
  #beta3 = splinefun(
  #  seq(0, 1, 0.125),
  #  c(0.025, 0.03, 0.04, 0.04, 0, -.17, -.35, -.47, -.52))

  # Yep, they look like the curves in the paper
  #xs = seq(0, 1, length.out = 101)
  #plot(xs, beta1(xs), type = "l", col = "blue", ylim = range(beta.true))
  #lines(xs, beta2(xs), col = "red")
  #lines(xs, beta3(xs), col = "green")

  # generate the X values
  X = array(rnorm(nT*nP, sd = sqrt(10)), c(nT, nP))
  tau = seq(0, 1, length.out = nD)

  # generate the mean value
  Y_mean = X[,1:3] %*% rbind(beta1(tau), beta2(tau), beta3(tau))

  # add errors
  Y = Y_mean + mvrnorm(nT, rep(0, nD), cov_mat) #rnorm_cork(nT, rep(0, nD), cov_mat)

  list(
    Y = Y,
    X = X,
    tau = tau,
    beta1 = beta1,
    beta2 = beta2,
    beta3 = beta3,
    Y_mean = Y_mean)
}

