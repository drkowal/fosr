library("fosr")

# sourcing files needed by iterate.R here
source("../R/helper_functions.R")
source("extra_mcmc_sampler.R")
source("modified_refund_lasso.R")
source("iterate.R")

Niters = 2

# Simulation parameters
n = 100     # Number of curves
m = 30      # Number of observation points
p_0 = 5    # Number of true zeros
p_1 = 5    # Number of true nonzeros (excluding the intercept)

# Correlation among predictors:
corr = 0.75

p = p_0 + p_1 + 1

gen_fosr_init = function()
{
  return(list(
    pci_cover_reg = array(dim = c(Niters, p)),
    pci_width_reg = array(dim = c(Niters, p)),
    MSE_reg = array(dim = c(Niters, p, 2), dimnames = list(1:Niters, 1:p, c("PM", "DSS"))),
    MSE_y = array(dim = c(Niters, 2), dimnames = list(1:Niters, c("PM", "DSS"))),
    TP = array(dim = c(Niters, p_0 + p_1, 2), dimnames = list(1:Niters, 1:(p_0+p_1), c("PM", "DSS"))),
    FP = array(dim = c(Niters, p_0 + p_1, 2), dimnames = list(1:Niters, 1:(p_0+p_1), c("PM", "DSS")))))
}

info_fosr       = gen_fosr_init()
info_fosr_basis = gen_fosr_init()
info_fosr_fpca  = gen_fosr_init()

info_refund_lasso = list(
  MSE_reg = array(dim = c(Niters, p)),
  MSE_y   = array(dim = Niters),
  TPR     = array(dim = Niters),
  TNR     = array(dim = Niters),
  TP      = array(dim = c(Niters, c(p_0 + p_1))),
  FP      = array(dim = c(Niters, c(p_0 + p_1))))

info_refund_ls = list(
  MSE_reg = array(dim = c(Niters, p)),
  MSE_y   = array(dim = Niters))

info_refund_gibbs = list(
  MSE_reg = array(dim = c(Niters, p)),
  MSE_y   = array(dim = Niters),
  pci_cover_reg = array(dim = c(Niters, p)),
  pci_width_reg = array(dim = c(Niters, p)))

for(ni in 1:Niters)
{
  set.seed(ni)

  sim_data = simulate_fosr(n = n, m = m, p_0 = p_0, p_1 = p_1,
                           RSNR = 5, K_true = 4,
                           sparse_factors = TRUE,
                           corr = corr,
                           perc_missing = 0)

  print("fosr")
  out = iterate(sim_data, "fosr")
  info_fosr$pci_cover_reg[ni,] = out$pci_cover_reg
  info_fosr$pci_width_reg[ni,] = out$pci_width_reg
  info_fosr$MSE_reg[ni,,] = out$MSE_reg
  info_fosr$MSE_y[ni,] = out$MSE_y
  info_fosr$TP[ni,,] = out$TP
  info_fosr$FP[ni,,] = out$FP

  print("fosr-basis")
  out = iterate(sim_data, "fosr-basis")
  info_fosr_basis$pci_cover_reg[ni,] = out$pci_cover_reg
  info_fosr_basis$pci_width_reg[ni,] = out$pci_width_reg
  info_fosr_basis$MSE_reg[ni,,] = out$MSE_reg
  info_fosr_basis$MSE_y[ni,] = out$MSE_y
  info_fosr_basis$TP[ni,,] = out$TP
  info_fosr_basis$FP[ni,,] = out$FP

  print("fosr-fpca")
  out = iterate(sim_data, "fosr-fpca")
  info_fosr_fpca$pci_cover_reg[ni,] = out$pci_cover_reg
  info_fosr_fpca$pci_width_reg[ni,] = out$pci_width_reg
  info_fosr_fpca$MSE_reg[ni,,] = out$MSE_reg
  info_fosr_fpca$MSE_y[ni,] = out$MSE_y
  info_fosr_fpca$TP[ni,,] = out$TP
  info_fosr_fpca$FP[ni,,] = out$FP

  print("refund:Lasso")
  out = iterate(sim_data, "refund:Lasso")
  info_refund_lasso$MSE_reg[ni,] = out$MSE_reg
  info_refund_lasso$MSE_y[ni] = out$MSE_y
  info_refund_lasso$TPR[ni] = out$TPR
  info_refund_lasso$TNR[ni] = out$TNR
  info_refund_lasso$TP[ni,] = out$TP
  info_refund_lasso$FP[ni,] = out$FP

  print("refund:GLS")
  out = iterate(sim_data, "refund:GLS")
  info_refund_ls$MSE_reg[ni,] = out$MSE_reg
  info_refund_ls$MSE_y[ni] = out$MSE_y

  print("refund:Gibbs")
  out = iterate(sim_data, "refund:Gibbs")
  info_refund_gibbs$MSE_reg[ni,] = out$MSE_reg
  info_refund_gibbs$MSE_y[ni] = out$MSE_y
  info_refund_gibbs$pci_cover_reg[ni,] = out$pci_cover_reg
  info_refund_gibbs$pci_width_reg[ni,] = out$pci_width_reg
}

#saveRDS(list(
#  fosr = info_fosr,
#  fosr_basis = info_fosr_basis,
#  fosr_fpca = info_fosr_fpca,
#  refund_lasso = info_refund_lasso,
#  refund_ls = info_refund_ls,
#  refund_gibbs = info_refund_gibbs), "simulation_out.rds")


