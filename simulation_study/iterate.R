source("../R/helper_functions.R")
source("extra_mcmc_sampler.R") # for fosr_basis
source("modified_refund_lasso.R") # for modified_fosr.vs_for_lasso
library("coda")
library("fosr")

iterate = function(
    sim_data, # the data to use
    type) # which method to use. One of
          # "fosr", "fosr-fpca", "fosr-basis",
          # "refund:Lasso", "refund:GLS", "refund:Gibbs"
{
  if(type == "fosr" || type == "fosr-fpca" || type == "fosr-basis")
  {
    # Coverage probabilities (indicators) and widths:
    pci_cover_reg = pci_width_reg = array(0, p)

    # Mean squared error:
    MSE_reg = array(0, c(p, 2), dimnames = list(1:p, c("PM", "DSS")))
    MSE_y = array(0, c(2), dimnames = list(c("PM", "DSS")))

    # True positive count and true negative count
    TP = array(NA, c(p_0 + p_1, 2)) # niter x model_size x c("DSS", "GBPV")
    FP = array(NA, c(p_0 + p_1, 2)) # niter x model_size x c("DSS", "GBPV")

    # Data:
    Y = sim_data$Y; X = sim_data$X; tau = sim_data$tau

    # Dimensions:
    n = nrow(Y); m = ncol(Y); p = ncol(X)
    #----------------------------------------------------------------------------
    # Run the FOSR:
    #----------------------------------------------------------------------------
    if(type == "fosr")
    {
      out = fosr::fosr(Y = Y, tau = tau, X = X,
                 K = 6, mcmc_params = list("beta", "fk", "alpha", "Yhat", "sigma_e", "sigma_g"))
    } else {
      use_fpca = (type == "fosr-fpca")
      out = fosr_basis(Y = Y, tau = tau, X = X,
                       use_fpca = use_fpca,
                       mcmc_params = list("beta", "fk", "alpha", "Yhat", "sigma_e", "sigma_g"))
    }
    #----------------------------------------------------------------------------
    # Decoupled shrinkage and selection:
    #----------------------------------------------------------------------------
    dss_info = fosr_select(X = X,
                            post_alpha = out$alpha,
                            post_trace_sigma_2 = n*m*out$sigma_e^2 + apply(out$sigma_g^2, 1, sum),
                            weighted = TRUE,
                            alpha_level = 0.10,
                            remove_int = TRUE,
                            include_plot = FALSE,
                            include_model_list = TRUE)
    alpha_dss = dss_info$alpha_dss
    model_list = dss_info$model_list[,-1] # remove intercept

    # remove intercept
    pos_select_true = which(
      apply(sim_data$alpha_tilde_true[,-1], 2, function(x) any(x != 0)))

    # Set TP and FP counts for DSS
    for(idx in 1:(dim(model_list)[1]))
    {
      model = which(model_list[idx,])
      model_size = length(model)

      TP[model_size,1] = sum(model %in% pos_select_true)
      FP[model_size,1] = model_size - TP[model_size,1]
    }

    # Set TP and FP counts for GBPV
    GBPV = fosr_gbpv(post_fk = out$fk, post_alpha = out$alpha)[-1] # remove intercept

    for(alpha in unique(GBPV))
    {
      model = which(GBPV <= alpha)
      model_size = length(model)

      TP[model_size,2] = sum(model %in% pos_select_true)
      FP[model_size,2] = model_size - TP[model_size,2]
    }

    #----------------------------------------------------------------------------
    # MSEs:

    # Yhat:
    MSE_y["PM"] = mean((sim_data$Y_true - colMeans(out$Yhat))^2)

    # Posterior mean:
    alpha_tilde_hat = matrix(0, nrow = m, ncol = p)
    Nsims = nrow(out$fk);
    for(nsi in 1:Nsims) {
      alpha_tilde_hat = alpha_tilde_hat + 1/Nsims*tcrossprod(out$fk[nsi,,], out$alpha[nsi,,])
    }

    # DSS estimate:
    alpha_tilde_dss = tcrossprod(colMeans(out$fk), alpha_dss)

    # DSS Yhat:
    MSE_y["DSS"] = mean(colMeans((sim_data$alpha_tilde_true - alpha_tilde_dss)^2))

    # MSEs:
    MSE_reg[,1] = colMeans((sim_data$alpha_tilde_true - alpha_tilde_hat)^2)
    MSE_reg[,2] = colMeans((sim_data$alpha_tilde_true - alpha_tilde_dss)^2)
    #----------------------------------------------------------------------------
    # Coverage probabilities:
    for(j in 1:p){
      post_alpha_tilde_j = get_post_alpha_tilde(post_fk = out$fk,
                                                post_alpha_j = out$alpha[,j,])
      ci_j = HPDinterval(as.mcmc(post_alpha_tilde_j))
      pci_cover_reg[j] = mean((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
      pci_width_reg[j] = mean(ci_j[,2] - ci_j[,1])
    }
    return(list(pci_cover_reg=pci_cover_reg,pci_width_reg=pci_width_reg,
             TP=TP,FP=FP,
                MSE_reg=MSE_reg,
                MSE_y=MSE_y))
  } else { # refund method
    # Coverage probabilities (indicators) and widths for gibbs
    pci_cover_reg = pci_width_reg = array(0, p)

    # Mean squared error:
    MSE_reg = array(0, p)
    MSE_y = array(0, 1)

    # True Positive Rate (sensitivity) and True Negative Rate (specificity):
    # for lasso
    TPR = TNR = array(0, 1)
    TP = array(NA, c(p_0 + p_1))
    FP = array(NA, c(p_0 + p_1))

    # Data:
    Y = sim_data$Y; X = sim_data$X; tau = sim_data$tau

    # Dimensions:
    n = nrow(Y); m = ncol(Y); p = ncol(X)
    data = as.data.frame(X); data$Y = Y;

    if(type == "refund:Lasso")
    {
      fit = refund::fosr.vs(Y ~ X - 1, data = data, method = "grLasso")
      pos_select_lasso = which(apply(t(fit$coefficients), 2, function(x) any(x != 0)))
      pos_select_true = which(apply(sim_data$alpha_tilde_true, 2, function(x) any(x != 0)))

      neg_select_lasso = (1:p)[-pos_select_lasso]
      neg_select_true  = (1:p)[-pos_select_true]

      MSE_y = mean((sim_data$Y_true - fit$fitted.values)^2)
      MSE_reg = colMeans((sim_data$alpha_tilde_true - t(fit$coefficients))^2)

      TPR = sum(pos_select_lasso %in% pos_select_true) / length(pos_select_true)
      TNR = sum(neg_select_lasso %in% neg_select_true) / length(neg_select_true)

      # Calculate ROC curve information now
      try_lambdas = exp(seq(log(1e-4), log(2), length.out = 150))
      try_lambdas_orig = try_lambdas
      msizes = c()
      while(length(try_lambdas) != 0)
      {
        ntry = length(try_lambdas)
        lambda.idx = ceiling(ntry/2)
        lambda = try_lambdas[lambda.idx]

        fit.info = modified_fosr.vs_for_lasso(
          Y ~ X - 1,
          data = data,
          lambda = lambda,
          method = "grLasso")

        model = which(apply(t(fit.info$coefficients)[,-1], 2, function(x) any(x != 0)))
        model_size = length(model)

        if(model_size == 0)
        {
          if(lambda.idx == 1)
            try_lambdas = c()
          else
            try_lambdas = try_lambdas[1:(lambda.idx-1)]
        }
        else if(model_size == p_0 + p_1)
        {
          if(lambda.idx == ntry)
            try_lambdas = c()
          else
            try_lambdas = try_lambdas[(lambda.idx+1):ntry]
        }
        else #record information
        {
          TP[model_size] = sum(model %in% pos_select_true)
          FP[model_size] = model_size - TP[model_size]

          # remove tihs lambda only
          try_lambdas = try_lambdas[-lambda.idx]

          # just want to make sure we get enough different model sizes
          msizes = c(msizes, model_size)
        }
      }
      return(list(MSE_y=MSE_y,MSE_reg=MSE_reg,TPR=TPR,TNR=TNR,TP=TP,FP=FP))
    } else if(type == "refund:GLS")
    {
      fit = refund::fosr.vs(Y ~ X - 1, data = data, method="ls")
      MSE_y = mean((sim_data$Y_true - fit$fitted.values)^2)
      MSE_reg = colMeans((sim_data$alpha_tilde_true - t(fit$coefficients))^2)
      return(list(MSE_y=MSE_y,MSE_reg=MSE_reg))
    } else if(type == "refund:Gibbs")
    {
      fit = refund::bayes_fosr(Y ~ X - 1, data = data, est.method = "Gibbs",
                               N.iter = 2000, N.burn = 1000)
      MSE_y = mean((sim_data$Y_true - fit$Yhat)^2)
      MSE_reg = colMeans((sim_data$alpha_tilde_true - t(fit$beta.hat))^2)
      #----------------------------------------------------------------------------
      # Coverage probabilities:
      for(j in 1:p){
        ci_j = cbind(fit$beta.LB[j,],
                     fit$beta.UB[j,])
        pci_cover_reg[j] = mean((sim_data$alpha_tilde_true[,j] >= ci_j[,1])*(sim_data$alpha_tilde_true[,j] <= ci_j[,2]))
        pci_width_reg[j] = mean(ci_j[,2] - ci_j[,1])
      }
      return(list(MSE_y=MSE_y,MSE_reg=MSE_reg,pci_cover_reg=pci_cover_reg,pci_width_reg=pci_width_reg))
    }
  }
}
