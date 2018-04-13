#source("dfosr_source/component_samplers.R") # for sampleMGP

# - If a1 is null, place a gamma(2, 1) prior on a1 otherwise the value
# - If a2 is null, place a gamma(2, 1) prior on a2 otherwise the value
# - h is the dimensions of delta to use
# - trans is the function that grabs theta.jh and delta.h for sampleMGP:
#     trans(this, up) returns list(theta.jh = ..., delta.h = ...)
gibber_mgp = function(
  trans,
  h,
  name,
  init = NULL,
  a1 = NULL,
  a2 = NULL)
{
  trans_mgp_a = function(this, up)
  {
    delta = up$val
    aa = this$val
    return(list(delta, aa))
  }

  sampleMGP_a1 = function(delta, a1_prev)
  {
    log_prob = function(a)
    {
      dgamma(delta[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)
    }

    uni.slice(a1_prev, g = log_prob, lower = 0, upper = Inf)
  }

  sampleMGP_a2 = function(delta, a2_prev)
  {
    log_prob = function(a)
    {
      sum(dgamma(delta[-1], shape = a, rate = 1, log = TRUE)) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)
    }

    uni.slice(a2_prev, g = log_prob, lower = 0, upper = Inf)
  }

  name_a1 = paste0("a1_", name)
  name_a2 = paste0("a2_", name)

  if(is.null(a1))
  {
    gib_a1 = gibber(
      sampleMGP_a1,
      2,
      list(),
      name_a1,
      trans_mgp_a,
      time = FALSE)
  }
  else
  {
    gib_a1 = gibber(function(a, b){ a1 }, a1, list(), name_a1, time = FALSE)
  }

  if(is.null(a1))
  {
    gib_a2 = gibber(
      sampleMGP_a2,
      3,
      list(),
      name_a2,
      trans_mgp_a,
      time = FALSE)
  }
  else
  {
    gib_a2 = gibber(function(a, b){ a2 }, a2, list(), name_a2, time = FALSE)
  }

  sample_this = gibber::compose_l(
    sampleMGP,
    function(this, up)
    {
      ll = trans(this, up)
      ll$a1 = get_val(this, name_a1)
      ll$a2 = get_val(this, name_a2)
      ll
    })

  if(is.null(init))
  {
    init = rep(1, h)
  }

  gibber(
    sample_this,
    init,
    list(gib_a1, gib_a2),
    name)
}

# omega is a n x p matrix
# omega_ij sim N(0, tau_ij^-1)
# tau_ij| x_ij sim G(1/2, x_ij)
# x_ij sim G(1/2, 1/sigma_e^2)
gibber_evol_independent_horseshoe = function(
  trans,
  n,
  p,
  name,
  init = NULL,
  sigma_e = 1/sqrt(n))
{
  if(is.null(init))
  {
    init = array(rgamma(n*p, 1, 1/sigma_e^2), dim = c(n, p))
  }

  name_x = paste0(name, "_x")

  gib_x = gibber(
    function(this, up)
    {
      tau = 1/(up$val^2)
      array(rgamma(n*p, 1, 1/sigma_e^2 + tau), dim = c(n, p))
    },
    init, # this value is just needed for the dimension
    list(),
    name_x,
    time = FALSE)

  sample_sig = function(omega, x)
  {
    tau = array(rgamma(n*p, 1, x + omega^2/2), dim = c(n, p))
    1/sqrt(tau)
  }

  sample_this = gibber::compose_l(
    sample_sig,
    function(this, up)
    {
      ll = trans(this, up)
      ll$x = get_val(this, name_x)
      ll
    })

  gibber(
    sample_this,
    init,
    list(gib_x),
    name)
}

# omega_np sim N(0, sig_p^2)
# sig_p^2 = 1/tau_p
# tau_p sim G(alpha, beta)
# p is the dimension of the gibber
gibber_normal_inverse_gamma = function(
  trans,
  p,
  name,
  init = NULL, # p x 1
  alpha = 0.01,
  beta = 0.01)
{
  # omega is a n x p matrix for any n
  #
  # omega_np sim N(0, tau_p^-1)
  # tau_p sim G(alpha, beta)
  sample_nig = function(omega)
  {
    n = nrow(omega)

    apply(
      omega,
      2,
      function(x)
      {
        post_shape = n/2 + alpha
        post_rate  = sum(x^2)/2 + beta
        1/sqrt(rgamma(1, shape = post_shape, rate = post_rate))
      })
  }

  if(is.null(init))
  {
    init = 1/sqrt(rgamma(p, max(0.01, alpha), max(0.01, beta)))
  }

  gibber(
    sample_nig,
    init,
    list(),
    name,
    trans)
}

#########################################################################

# omega is a n x p matrix
# omega_ij sim N(0, sig_ij^2)
# sig_ij = sr_i*sc_j
# sr_i^2 sim inverse gamma
# sc_j^2 = cumprod(deltac)
# deltac sim MGP
#
# trans maps (this, up) to omega
gibber_ig_row_mgp_col = function(
  trans,
  n,
  p,
  name,
  init,
  sigma_e = 1/sqrt(n))
{
  # want to return gibber that samples from s2u_ij given omega

  name_sr = paste0(name, "_sr")
  name_deltac = paste0(name, "_deltac")

  trans_sr = function(this, up)
  {
    omega = up$pass_val[[1]]
    sc = 1/sqrt(cumprod(get_val(up, name_deltac)))

    list(t(omega)*1/sc)
  }

  gib_sr = gibber_normal_inverse_gamma(
    trans_sr,
    n,
    name_sr)

  trans_deltac = function(this, up)
  {
    omega = up$pass_val[[1]]
    sr = get_val(up, name_sr)
    deltac = this$val

    list(omega*1/sr, deltac)
  }

  gib_deltac = gibber_mgp(
    trans_deltac,
    p,
    name_deltac)

  gibber(
    function(this, up)
    {
      sr = get_val(this, name_sr)
      sc = 1/sqrt(cumprod(get_val(this, name_deltac)))

      sig = tcrossprod(sr, sc)

      is_bad = which(
        sig < sqrt(.Machine$double.eps),
        arr.ind = TRUE)
      sig[is_bad] = sqrt(.Machine$double.eps)

      return(sig)
    },
    init,
    down = list(gib_sr, gib_deltac),
    name = name,
    time = FALSE,
    pass = trans) # omega gets set in pass_val before sr and sc are called
}

# omega is a n x p matrix
# omega_ij sim N(0, sig_ij^2)
# sig_ij sim C+(0, lambda_j)
# lambda_j sim C+(0, lambda_o)
# lambda_o sim C+(0, sig_e)
#
# The parameter expression goes like this
#  eta2 -> gam2 -> eta1 -> gam1 -> eta0 -> gam0
#  etaj is j dimensions
#
# trans maps (this, up) to omega
gibber_evol_column_horseshoe = function(
  trans,
  nN,
  nP,
  name,
  init,
  sigma_e = 1/sqrt(nN))
{
  sample_gam0 = function(this, up)
  {
    eta0 = up$val

    rgamma(1, 1, eta0 + 1/sigma_e^2)
  }

  gib_gam0 = gibber(
    sample_gam0,
    1,
    list(),
    "gam0",
    time = FALSE)

  sample_eta0 = function(this, up)
  {
    gam1 = up$val
    gam0 = get_val(this, "gam0")

    rgamma(1, (nN + 1)/2, gam0 + sum(gam1))
  }

  gib_eta0 = gibber(
    sample_eta0,
    1,
    list(gib_gam0),
    "eta0",
    time = FALSE)

  sample_gam1 = function(this, up)
  {
    eta1 = up$val
    eta0 = get_val(this, "eta0")

    rgamma(nN, 1, eta1 + eta0)
  }

  gib_gam1 = gibber(
    sample_gam1,
    rep(1, nN),
    list(gib_eta0),
    "gam1",
    time = FALSE)

  sample_eta1 = function(this, up)
  {
    gam2 = up$val
    gam1 = get_val(this, "gam1")

    rgamma(nN, (nP+1)/2, gam1 + rowSums(gam2))
  }

  gib_eta1 = gibber(
    sample_eta1,
    rep(1, nN),
    list(gib_gam1),
    "eta1",
    time = FALSE)

  sample_gam2 = function(this, up)
  {
    eta2 = (1/up$val)^2
    eta1 = get_val(this, "eta1")

    gam2 = array(dim = c(nN, nP))
    for(nidx in 1:nN)
    {
      gam2[nidx,] = rgamma(nP, rep(1, nP), eta2[nidx,] + eta1[nidx])
    }

    return(gam2)
  }

  gib_gam2 = gibber(
    sample_gam2,
    array(1, dim = c(nN, nP)),
    list(gib_eta1),
    "gam2",
    time = FALSE)

  trans_sig = function(this, up)
  {
    ll = trans(this, up) # this gives omega
    ll$gam2 = get_val(this, "gam2")
    return(ll)
  }

  sample_sig = function(
    omega,  # n x p
    gam2)   # n x p
  {
    eta = array(dim = c(nN, nP))
    for(p in 1:nP)
    {
      eta[,p] = rgamma(1, rep(1, nN), omega[,p]^2/2 + gam2[,p])
    }

    return(1/sqrt(eta))
  }

  gibber(
    sample_sig,
    init,
    list(gib_gam2),
    name,
    trans_sig,
    time = TRUE)
}

# omega is a n x p matrix
# omega_ij sim N(0, sig_ij^2)
# sig_ij sim C+(0, lambda_i)
# lambda_i sim C+(0, lambda_o)
# lambda_o sim C+(0, sig_e)
#
# The parameter expression goes like this
#  eta2 -> gam2 -> eta1 -> gam1 -> eta0 -> gam0
#  etaj is j dimensions
#
# trans maps (this, up) to omega
gibber_evol_row_horseshoe = function(
  trans,
  nN,
  nP,
  name,
  init,
  sigma_e = 1/sqrt(nN))
{
  name_gib_col = paste0(name, "_transpose")
  gib_col = gibber_evol_column_horseshoe(
    function(this, up)
    {
      up$pass_val # expose t(omega)
    },
    nP,
    nN,
    name_gib_col,
    t(init),
    sigma_e)

  gibber(
    function(this, up) # sample function
    {
      t(get_val(this, name_gib_col))
    },
    init,
    list(gib_col),
    name,
    time = FALSE,
    pass = function(this, up)
    {
      list(t(trans(this, up)[[1]]))  # get omega and transpose it
                                     # and this will be passed to gib_col
    })
}

