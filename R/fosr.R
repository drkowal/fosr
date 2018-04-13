# Y[t,] = Fmat%*%Beta[t,] + epsilon[t,]
# Beta[,k] = mu[,k] + X%*%alpha[,k] + gamma_tk[,k]
fosr = function(
  Y,
  tau,
  X = NULL,
  K = NULL,
  nsave = 1000,
  nburn = 1000,
  nskip = 3,
  mcmc_params = c("Beta", "Y_hat", "Fmat", "alpha_pk"),
  f_gibber_sig_alpha_pk = gibber_evol_column_horseshoe,
  f_gibber_sig_gamma_tk = gibber_evol_col_row)
{
  start_time = Sys.time()

  #############################################################################
  # -  define constants used throughout + initial values
  # -  define gibbers
  #      transform function,
  #      sample function,
  #      initial value
  # -  build model from gibber components
  # -  run gibbs sampler and return the output
  #############################################################################

  # not handling missing values at the moment
  stopifnot(!any(is.na(Y)))
  stopifnot(!any(is.na(tau)))

  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)

  # Compute the dimensions:
  nT = nrow(Y)
  nM = ncol(Y)
  nD = ncol(tau)
  nNotNA = sum(!is.na(Y))

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  fdlm_inits = fdlm_init_d(Y, tau, K)
  init_Beta  = fdlm_inits$Beta
  init_Psi   = fdlm_inits$Psi
  splineInfo = fdlm_inits$splineInfo

  Bmat  = splineInfo$Bmat
  Omega = splineInfo$Omega
  nK = ncol(init_Beta) # to be sure we have the right value

  BtY = tcrossprod(t(Bmat), Y)

  # FLC matrix:
  init_Fmat = Bmat%*%init_Psi

  # Initialize the conditional expectation:
  init_BetaPsit = tcrossprod(init_Beta, init_Psi)
  init_Btheta = tcrossprod(init_BetaPsit, Bmat)

  # Initialize the observation error SD:
  init_sig_eps = sd(Y - init_Btheta, na.rm=TRUE)
  init_sig_eps_t = rep(init_sig_eps, nT)

  # Initialize the FLC smoothing parameters (conditional MLE):
  init_tau_f_k = apply(
    init_Psi,
    2,
    function(x)
    {
      (ncol(Bmat) - (nD+1)) /
        crossprod(x, Omega)%*%x
    })

  # Predictors:
  if(!is.null(X))
  {
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
  X = cbind(rep(1, nT), X)

  # Number of predictors (including the intercept)
  nP = ncol(X)

  # Initialize the regression terms (and the mean term)
  init_alpha_pk = matrix(0, nrow = nP, ncol = nK) # Regression coefficients
  init_gamma_tk = matrix(0, nrow = nT, ncol = nK) # Residuals

  # Initialize the regression coefficients via sampling (p >= T) or OLS (p < T)
  for(k in 1:nK)
  {
    if(nP >= nT)
    {
      init_alpha_pk[,k] = sampleFastGaussian(
        Phi = X/init_sig_eps_t,
        Ddiag = rep(.01*init_sig_eps^2, nP),
        alpha = tcrossprod(Y, t(init_Fmat[,k]))/init_sig_eps)
    }
    else
    {
      init_alpha_pk[,k] = lm(init_Beta[,k] ~ X - 1)$coef
    }

    # Residuals:
    init_gamma_tk[,k] = init_Beta[,k] - X%*%init_alpha_pk[,k]
  }

  # Intercept term:
  init_mu_k = as.matrix(init_alpha_pk[1,])

  # SD term for mu_k:
  init_a1_mu = 2; init_a2_mu = 3
  init_delta_mu_k = sampleMGP(
    matrix(init_mu_k, ncol = nK),
    rep(1, nK),
    a1 = init_a1_mu,
    a2 = init_a2_mu)
  init_sig_mu_k = 1/sqrt(cumprod(init_delta_mu_k))

  # Initialize the corresponding SD term(s):
  evolParams_int = # this isn't used except for initial values
    initEvolParams(omega = init_gamma_tk, evol_error = "NIG")

  # Update the error SD for gamma:
  init_sig_gamma_t = evolParams_int$sigma_wt[,1]
  init_sig_gamma_tk = evolParams_int$sigma_wt*rep(rep(1, nK), each = nT)

  if(nP > 1)
  {
    evolParams = initEvolParams(
      omega = init_alpha_pk[-1,],
      evol_error = "HS")

    # SD for omega:
    init_sig_alpha_pk = evolParams$sigma_wt
  }

  #############################################################################
  # define gibbers
  #############################################################################

  #############################################################################
  # delta_mu_k

  # trans needs to define a list with the following values:
  #   theta.jh, delta.h
  # trans: this, up -> theta.jh, delta.h
  trans_delta_mu_k = function(this, up)
  {
    alpha_pk = up$val
    mu_k = alpha_pk[1,]

    list(
      theta.jh = matrix(mu_k, ncol = nK),
      delta.h = this$val)
  }

  gib_delta_mu_k = gibber_mgp(
    trans_delta_mu_k,
    nK,
    "delta_mu_k",
    init_delta_mu_k)

  #############################################################################
  # sig_mu_k
  sample_sig_mu_k = function(this, up)
  {
    delta_mu_k = get_val(up, "delta_mu_k");
    1/sqrt(cumprod(delta_mu_k))
  }

  gib_sig_mu_k = gibber(
    sample_sig_mu_k,
    init_sig_mu_k,
    list(),
    "sig_mu_k",
    time = FALSE)

  #############################################################################
  # sig_alpha_pk
  if(nP == 1)
  {
    gib_sig_alpha_pk = gibber_constant(NULL, list(), "sig_alpha_pk")
  }
  else
  {
    trans_sig_alpha_pk = function(this, up)
    {
      alpha = up$val
      list(array(alpha[-1,], c(nP-1, nK)))
    }

    gib_sig_alpha_pk = f_gibber_sig_alpha_pk(
      trans_sig_alpha_pk,
      nP-1,
      nK,
      "sig_alpha_pk",
      init_sig_alpha_pk)
  }

  #############################################################################
  # alpha_pk
  sample_alpha_pk = function(
    Y_tilde,
    sig_tilde,
    sig_gamma_tk,
    sig_alpha_pk,
    sig_mu_k)
  {
    alpha_pk = array(dim = c(nP, nK))

    # Draw Separately for each k:
    for(k in 1:nK)
    {
      # Marginalize over gamma_{tk} to sample {alpha_pk}_p for fixed k:
      y_tilde_k = Y_tilde[,k];
      sig_tilde_k = sqrt(sig_tilde^2 + sig_gamma_tk[,k]^2)

      if(nP >= nT)
      {
        # Fast sampler for p >= T (BHATTACHARYA et al., 2016)
        alpha_pk[,k] = sampleFastGaussian(
          Phi = X/sig_tilde_k,
          Ddiag = as.numeric(c(sig_mu_k[k], sig_alpha_pk[,k])^2),
          alpha = y_tilde_k/sig_tilde_k)
      }
      else
      {
        chQ_k = chol(
          crossprod(X/sig_tilde_k) +
          diag(as.numeric(1/c(sig_mu_k[k], sig_alpha_pk[,k])^2)))
        ell_k = crossprod(X, y_tilde_k/sig_tilde_k^2)

        alpha_pk[,k] = backsolve(
          chQ_k,
          forwardsolve(t(chQ_k), ell_k) + rnorm(nP))
      }
    }

    return(alpha_pk)
  }

  trans_alpha_pk = function(this, up)
  {
    list(
      Y_tilde      = get_val(up,   "Y_tilde"),
      sig_tilde    = get_val(up,   "sig_eps_t"),
      sig_gamma_tk = get_val(up, c("gamma_tk", "sig_gamma_tk")),
      sig_alpha_pk = get_val(this, "sig_alpha_pk"),
      sig_mu_k     = get_val(this, "sig_mu_k"))
  }

  gib_alpha_pk = gibber(
    sample_alpha_pk,
    init_alpha_pk,
    list(
      gib_delta_mu_k,
      gib_sig_mu_k,
      gib_sig_alpha_pk),
    "alpha_pk",
    trans_alpha_pk)

  #############################################################################
  # sig_gamma_t
  trans_sig_gamma_tk = function(this, up)
  {
    gamma_tk = up$val
    list(gamma_tk)
  }

  gib_sig_gamma_tk = f_gibber_sig_gamma_tk(
    trans_sig_gamma_tk,
    nT,
    nK,
    "sig_gamma_tk",
    init_sig_gamma_tk)

  #############################################################################
  # gamma_tk
  sample_gamma_tk = function(this, up)
  {
    Y_tilde      = get_val(up,   "Y_tilde")
    alpha_pk     = get_val(up,   "alpha_pk")
    sig_tilde    = get_val(up,   "sig_eps_t")
    sig_gamma_tk = get_val(this, "sig_gamma_tk")

    post_sd   =
      1/sqrt(rep(1/sig_tilde^2, times = nK) + matrix(1/sig_gamma_tk^2))
    post_mean =
      matrix((Y_tilde - X%*%alpha_pk)/rep(sig_tilde^2, times = nK))*post_sd^2
    gamma_tk = array(
      rnorm(n = nT*nK, mean = post_mean, sd = post_sd),
      c(nT, nK))

    return(gamma_tk)
  }

  gib_gamma_tk = gibber(
    sample_gamma_tk,
    init_gamma_tk,
    list(gib_sig_gamma_tk),
    "gamma_tk")

  #############################################################################
  # sig_eps_t
  sample_sig_eps_t = function(this, up)
  {
    Y_hat = get_val(up, "Y_hat")

    post_shape = 0.5*nNotNA
    post_rate  = 0.5*sum((Y - Y_hat)^2, na.rm = TRUE)

    sig_eps = 1/sqrt(rgamma(1, shape = post_shape, rate = post_rate))

    rep(sig_eps, nT)
  }

  gib_sig_eps_t = gibber(
    sample_sig_eps_t,
    init_sig_eps_t,
    list(),
    "sig_eps_t")

  #############################################################################
  # Beta
  gib_Beta = gibber(
    function(this, up)
    {
      alpha_pk = get_val(up, "alpha_pk")
      gamma_tk = get_val(up, "gamma_tk")

      X %*% alpha_pk + gamma_tk
    },
    init_Beta,
    list(),
    "Beta",
    time = FALSE)

  #############################################################################
  # Y_hat
  gib_Y_hat = gibber(
    function(this, up)
    {
      Beta = get_val(up, "Beta")
      Fmat = get_val(up, "Fmat")

      tcrossprod(Beta, Fmat)
    },
    Y,
    list(),
    "Y_hat",
    time = FALSE)

  #############################################################################
  # Y_tilde
  gib_Y_tilde = gibber(
    function(this, up)
    {
      Psi = get_val(up, "Psi")
      crossprod(BtY, Psi)
    },
    Y %*% init_Fmat,
    list(),
    "Y_tilde",
    time = FALSE)

  #############################################################################
  # Fmat
  gib_Fmat = gibber(
    function(this, up)
    {
      Psi = get_val(up, "Psi")
      Bmat %*% Psi
    },
    init_Fmat,
    list(),
    "Fmat",
    time = FALSE)

  #############################################################################
  # sample tau_f_k
  trans_tau_f_k = function(this, up)
  {
    list(
      this$val,
      up$val,
      Omega = Omega,
      d = nD,
      uniformPrior = TRUE,
      orderLambdas = TRUE)
  }

  gib_tau_f_k = gibber(
    sample_lambda,
    init_tau_f_k,
    list(),
    "tau_f_k",
    trans_tau_f_k)

  #############################################################################
  # Psi
  trans_Psi = function(this, up)
  {
    list(
      BtY = BtY,
      Beta = get_val(up, "Beta"),
      Psi = this$val,
      BtB = diag(nrow(BtY)),
      Omega = Omega,
      lambda = get_val(this, "tau_f_k"),
      sigmat2 = get_val(up, "sig_eps_t")^2)
  }

  gib_Psi = gibber(
    fdlm_flc,
    init_Psi,
    list(gib_tau_f_k),
    "Psi",
    trans_Psi)

  #############################################################################
  # the model
  model = gibber_constant(
    NULL, # really, Y and the list of constants
          # that were defined earlier in this function
          # should be incorporated here
    list(
      gib_Psi,
      gib_Fmat,
      gib_Beta,
      gib_Y_hat,
      gib_Y_tilde,
      gib_alpha_pk,
      gib_sig_eps_t,
      gib_gamma_tk),
    "fosr")

  # run the gibbs sampler of the constructed model
  model = gibbs(
    model,
    keep = c(),
    nsave = 0,
    nburn = 10,
    nskip = 0)$gibber

  # adjust the ordering
  if(nK > 1)
  {
    tau_f_k = get_val(model, c("Psi", "tau_f_k"))
    Psi     = get_val(model, "Psi")
    Beta    = get_val(model, "Beta")

    adjOrder = order(tau_f_k, decreasing = TRUE);
    tau_f_k = tau_f_k[adjOrder];
    Psi = Psi[,adjOrder];
    Beta = as.matrix(Beta[,adjOrder])

    # The gibber construction is designed to make
    # modifications in place like this difficult
    # and not impossible.
    model$down[["Psi"]]$down[["tau_f_k"]]$val = tau_f_k
    model$down[["Psi"]]$val = Psi
    model$down[["Beta"]]$val = Beta
  }

  out = gibbs(
    model,
    keep  = mcmc_params,
    nsave = nsave,
    nburn = max(0, nburn - 10),
    nskip = nskip)

  end_time = Sys.time()
  out$run_time = difftime(end_time, start_time, units = "secs")
  out$nsave = nsave
  out$nburn = nburn
  out$nskip = nskip

  out
}
