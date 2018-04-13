
# mat: a matrix with values of a function to plot
# tau: the x values. If not given, assumed to be uniform on (0, 1)
# along: which dimension of mat contains the individual functions
plot_spaghetti = function(mat, tau = NULL, along = 1, main = "", xlab = "tau", ylab = "")
{
  nD = ifelse(along == 1, ncol(mat), nrow(mat))
  nF = ifelse(along == 1, nrow(mat), ncol(mat))

  if(is.null(tau))
    tau = seq(0, 1, length.out = nD)

  if(along == 1)
  {
    plot(
      tau, mat[1,], type = "l",
      ylim = range(mat),
      xlab = xlab, ylab = ylab,
      main = main)
  }
  else
  {
    plot(
      tau, mat[,1], type = "l",
      ylim = range(mat),
      xlab = xlab, ylab = ylab,
      main = main)
  }

  if(nF > 1)
  {
    for(idx in 2:nF)
    {
      if(along == 1)
        lines(tau, mat[idx,])
      else
        lines(tau, mat[,idx])
    }
  }
}

# given alpha and F, F*alpha gives a M x P matrix where
# M is the number of evaluation points and P is the number
# of predictors. So F*alpha is the P predictor functions
# in the model Y_t(s) = sum_j Xij predictor_j(s) + error_t(s)
#
# alpha: nMCMC x nP x nK matrix, from nMCMC draws from some density
# Fmat:  nMCMC x nM x nK matrix, from nMCMC draws from some density
#
# The return value is a nMCMC x nM x nP matrix
calc_post_predictor = function(alpha, Fmat)
{
  nMCMC = length(Fmat[,1,1])
  nM = nrow(Fmat[1,,])
  nP = nrow(alpha[1,,])

  post_g = array(dim = c(nMCMC, nM, nP))

  for(idx in 1:nMCMC)
  {
    post_g[idx,,] = tcrossprod(Fmat[idx,,], alpha[idx,,])
  }

  post_g
}

calc_mat_mse = function(draws, true)
{
  draws_dim = dim(draws)
  tmp = array(dim = draws_dim)
  nMCMC = draws_dim[1]

  for(idx in 1:nMCMC)
    tmp[idx,,] = (draws[idx,,] - true)^2

  meanAlong(tmp, 1)
}

# draws: nMCMC x nR x nC
# true:  nR x nC
#
# returns (draws[idx,,] - true)^2 for 1:nMCMC
calc_mat_mse_many = function(draws, true)
{
  dim_draws = dim(draws)
  ret = array(dim = dim_draws)

  for(idx in 1:dim_draws[1])
    ret[idx,,] = (draws[idx,,] - true)^2

  ret
}


