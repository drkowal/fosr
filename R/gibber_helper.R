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
# omega_ij        sim N(0, tau_ij^-1)
# tau_ij | x_ij   sim G(1/2, x_ij)
# x_ij | tt_j     sim G(1/2, tt_j)
# tt_j | xx_j     sim G(1/2, xx_j)
# xx_j            sim G(1/2, 1/sigma_e^2)
#
# trans maps (this, up) to omega
gibber_evol_column_horseshoe = function(
  trans,
  n,
  p,
  name,
  init = NULL,
  sigma_e = 1/sqrt(n))
{
  name_xx = paste0(name, "_xx")
  name_x  = paste0(name, "_x")
  name_tt = paste0(name, "_tt")

  gib_xx = gibber(
    function(this, up)
    {
      tt = up$val
      rgamma(p, 1, 1/sigma_e^2 + tt)
    },
    rgamma(p, 1, 1/sigma_e^2),
    list(),
    name_xx,
    time = FALSE)

  gib_tt = gibber(
    function(this, up)
    {
      xx = get_val(this, name_xx)
      x  = up$val

      rgamma(p, 0.5 + n/2, colSums(x) + xx)
    },
    rgamma(p, 0.5 + n/2, 1),
    list(gib_xx),
    name_tt,
    time = FALSE)

  gib_x = gibber(
    function(this, up)
    {
      tau = 1/(up$val^2)
      tt  = get_val(this, name_tt)

      array(
        rgamma(n*p, 1, tau + tcrossprod(rep(1, n), tt)),
        dim = c(n, p))
    },
    array(rgamma(n*p, 1, 1), dim = c(n, p)),
    list(gib_tt),
    name_x,
    time = FALSE)

  sample_sig = function(omega, x)
  {
    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(
      rep(1,n),
      apply(
        omega,
        2,
        function(ccc)
        {
          any(ccc^2 < 10^-16)*max(10^-8, mad(ccc)/10^6)
        }))
    hsInput2 = omega^2 + hsOffset

    tau = array(rgamma(n*p, 1, x + hsInput2/2), c(n, p))
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
gibber_evol_col_row = function(
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
