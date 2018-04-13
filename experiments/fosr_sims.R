# Libraries we need:
library(truncdist);
library(coda);
library(dsp);
library(BayesLogit);
library(fields);
library(Rcpp);
library(KFAS);
library(refund)
library(fds);

# Source files:
library("fosr", lib.loc = "~/ps/bhdfosr/fosr/bld")
##source('dfosr_source/commonSV_source.R')
##source('dfosr_source/component_samplers.R')
##source('dfosr_source/helper_functions.R')
##source('dfosr_source/mcmc_sampler.R')
##source('dfosr_source/extra_samplers.R')
##source('dfosr_source/extra_samplers_basis.R')
##Rcpp::sourceCpp('dfosr_source/fdlm_source.cpp')

#set.seed(14850)
#----------------------------------------------------------------------------
# Simulate some data:
#----------------------------------------------------------------------------
T = 200;   # Number of time points (or subjects)
m = 100;   # Number of observation points
RSNR = 5;  # Root signal-to-noise ratio
ar1 = 0;   # Autoregressive coefficient of the predictors X_j (i.e., autocorrelated?)
K_true = 4 # Number of factors
p_0 = 2;   # Number of zero coefficients
p_1 = 5;   # Number of nonzero "coefficients" (predictors j for which \alpha_{k,j} != 0 for some k)

# Number of predictors:
p = 1 + p_1 + p_0

# Observation points:
tau = seq(0, 1, length.out = m)

# FLCs, orthogonalized, are polynomials centered at 1/2
F_true = matrix(0, nrow = m, ncol = K_true)
for(k in 1:K_true) F_true[,k] = (tau - 0.5)^(k-1)
# Orthonormalize:
F_true = qr.Q(qr(F_true))
# Positive:
for(k in 1:K_true){if(sum(F_true[,k]) < 0) F_true[,k] = -F_true[,k]}
#plot(as.ts(F_true))

# Simulate the predictors:
if(ar1 == 0){
  X = cbind(1,matrix(rnorm(n = T*(p-1)), nr = T, nc = p-1))
} else X = cbind(1,
                 apply(matrix(0, nr = T, nc = p-1), 2, function(x)
                   arima.sim(n = T, list(ar = ar1), sd = sqrt(1-ar1^2))))
#plot(as.ts(X))

# True coefficients:
alpha.arr.true = array(0, c(p, K_true))

# p = 1 Intercept coefficients: decaying importance
alpha.arr.true[1,] = 1/(1:K_true)

# Randomly select which k's are nonzero (k_p) for each predictor:
if(p_1 > 0){for(j in 1:p_1){
  # Which factors are nonzero for predictor j? Truncated Poisson(1)
  k_p = sample(1:K_true, rtrunc(n = 1, spec = 'pois', b = K_true, lambda = 1))

  # Simulate true values of the (nonzero) regression coefficients (decaying importance)
  alpha.arr.true[j+1, k_p] = rnorm(n = length(k_p), mean = 0, sd = 1/k_p)
}}

# True factors: add Gaussian (subject-specific) noise (decaying importance)
Beta_true = matrix(0, nrow = T, ncol = K_true)
for(k in 1:K_true) Beta_true[,k] = X%*%alpha.arr.true[,k] + rnorm(n = T, sd = 1/k)

# True FTS:
Y_true = tcrossprod(Beta_true, F_true)

# Noise SD based on RSNR:
sigma_true = sd(Y_true)/RSNR

# Observed data:
Y = Y_true + sigma_true*rnorm(m*T)

# Plot the true curves:
#plot(fts(x = tau, y = t(Y_true)), plot.type = "functions")

# Plot the observed curves:
# plot(fts(x = tau, y = t(Y)), plot.type = "functions")

#----------------------------------------------------------------------------
# FOSR:
sss = Sys.time()
out = fosr(Y = Y, tau = tau, X = X,
           K = 6, dist_reg_coef = "HS", add_p_shrink = FALSE,
           mcmc_params = list("beta", "fk", "alpha", "Yhat"),
           nsave = 1000, nburn = 0, nskip = 0)
eee = Sys.time() # (5.08 w/ gibber.. no F.., 8.10 here@1000)
                 # 76@10000
difftime(eee, sss, units = "secs")

# Plot the fitted alpha's w/ HPD intervals:
alpha_hat = colMeans(out$alpha);  K = ncol(alpha_hat)
for(j in p:1){
  dev.new()
  ci_j = HPDinterval(as.mcmc(out$alpha[,j,]))
  plot(1:K, alpha_hat[j,], ylim = range(ci_j, alpha.arr.true[j,]), pch = paste(1:K), main = paste('Predictor j =',j), xlab = 'Factor Number', ylab = expression(alpha[kj]))
  for(k in 1:K) {lines(rep(k, 2), ci_j[k,], col='blue', lwd=4); if(k <= K_true) lines(k, alpha.arr.true[j,k], type='p', pch=4, lwd=2, cex = 2) }
}

# Compare the fitted regression functions:
g_true = tcrossprod(F_true, alpha.arr.true) # m x p
Ni = nrow(out$fk);
post_g = array(0, c(Ni, m, p));
for(ni in 1:Ni)
  post_g[ni,,] = tcrossprod(out$fk[ni,,], out$alpha[ni,,])
g_hat = colMeans(post_g)
for(j in p:1)
{
  dev.new()
  plot_curve( # helper_functions.R
    post_f = post_g[,,j],
    tau = tau,
    include_joint = TRUE,
    main = paste('Predictor',j));
  lines(tau, g_true[,j], lwd = 8, lty =3)
}

# Include some other plots:
plot_factors(out$beta)
plot_flc(out$fk, tau = tau)
i = sample(1:T, 1); plot_fitted(y = Y[i,], mu = colMeans(out$Yhat[,i,]), postY = out$Yhat[,i,], y_true = Y_true[i,],t01 = tau)

# Mean squared error:
(MSE0 = colMeans((g_true - g_hat)^2));
(RMSE0 = sqrt(mean(MSE0)))

# Existing Methods:
if(FALSE){
  out0 = naive_lm(Y, X)
  g_hat_naive = t(colMeans(out0$alpha))
  (MSE1 = colMeans((g_true - g_hat_naive)^2));
  (RMSE1 = sqrt(mean(MSE1)))

  # Compare these posterior distributions:
  for(j in p:1){
    plot_curve(post_f = post_g[,,j], tau = tau, include_joint = TRUE, main = paste('FOSR, Predictor',j)); lines(tau, g_true[,j], lwd = 8, lty =3)
    plot_curve(post_f = out0$alpha[,j,], tau = tau, include_joint = TRUE, main = paste('Naive, Predictor',j)); lines(tau, g_true[,j], lwd = 8, lty =3)
  }

  timer0 = proc.time()[3] # For timing the procedure
  data = as.data.frame(X); data$Y = Y
  #fit_fosr = refund::fosr(Y~X-1, data = data, Y = Y, X = X) #plot(fit_fosr_vs)
  fit_fosr_lasso = refund::fosr.vs(Y ~ X - 1, data = data, method = "grLasso") #plot(fit_fosr_vs)
  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))
  g_hat_lasso = t(fit_fosr_lasso$coefficients)

  (MSE2 = colMeans((g_true - g_hat_lasso)^2));
  (RMSE2 = sqrt(mean(MSE2)))

  fit_fosr_ls = refund::fosr.vs(Y ~ X - 1, data = data, method="ls") #plot(fit_fosr_vs)
  g_hat_ls = t(fit_fosr_ls$coefficients)


  (MSE3 = colMeans((g_true - g_hat_ls)^2));
  (RMSE3 = sqrt(mean(MSE3)))


  RMSE0; RMSE1; RMSE2; RMSE3
  #percentRedMSE = 100*(MSE1 - MSE0)/MSE1
  #plot(1:p, percentRedMSE, ylim = range(c(0, 100, percentRedMSE)), main = 'Percent reduction in MSE')
}

# Global Bayesian P-values:
(GBPV = apply(post_g, 3, function(x) min(simBaS(x))))

# ESS:
getEffSize(out$fk)
getEffSize(out$beta)
getEffSize(out$alpha)
getEffSize(post_g)

