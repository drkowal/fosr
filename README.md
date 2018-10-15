# fosr
Bayesian Function-on-Scalar Regression for High Dimensional Data

# Instructions
- Clone the code to `/path/to/fosr`
- `R CMD INSTALL /path/to/fosr`
- In R, load the library: `library(fosr)`

# Example usage

```
library(fosr)

# Simulate some data:
n = 100
m = 20
p_0 = 100
p_1 = 5
sim_data = simulate_fosr(n = n, m = m, p_0 = p_0, p_1 = p_1)

# Data:
Y = sim_data$Y
X = sim_data$X
tau = sim_data$tau

# Dimensions:
n = nrow(Y)
m = ncol(Y)
p = ncol(X)

# Run the FOSR:
out = fosr(
  Y = Y,
  tau = tau,
  X = X,
  K = 6,
  mcmc_params = list("fk", "alpha", "Yhat", "sigma_e", "sigma_g"))

# Plot a posterior summary of a regression function, say j = 3:
j = 3
post_alpha_tilde_j = get_post_alpha_tilde(out$fk, out$alpha[,j,])
plot_curve(post_alpha_tilde_j, tau = tau)

# Add the true curve:
lines(tau, sim_data$alpha_tilde_true[,j], lwd=6, col='green', lty=6)

# Plot the loading curves:
plot_flc(out$fk, tau = tau)

# Plot the fitted values for a random subject:
i = sample(1:n, 1)
plot_fitted(y = Y[i,], mu = colMeans(out$Yhat[,i,]),
            postY = out$Yhat[,i,], y_true = sim_data$Y_true[i,], t01 = tau)

# perform DSS variable selection
alpha_dss = fosr_select(
  X = X,
  post_alpha = out$alpha,
  post_trace_sigma_2 = n*m*out$sigma_e^2 + apply(out$sigma_g^2, 1, sum),
  weighted = TRUE,
  alpha_level = 0.10,
  remove_int = TRUE,
  include_plot = TRUE)

# Which variables do we include?
pos_select_dss = which(apply(alpha_dss, 1, function(x) any(x != 0)))

# And the truth
pos_select_true = which(apply(sim_data$alpha_tilde_true, 2, function(x) any(x != 0)))
```
