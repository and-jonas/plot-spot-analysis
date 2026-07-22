
# get full correlation matrix
get_cor_matrix <- function(corStruct, n){
  phi <- coef(corStruct, unconstrained = FALSE)
  mat <- phi ^ (abs(outer(1:n,1:n, "-")))
  return(mat)
}

# invert to original scale
inv_logit_adjusted <- function(x, eps=0.001){
  p <- exp(x)/(1+exp(x))
  p*(1+eps) - eps
}

# estimate effective sample size
n_eff_from_acf <- function(n, rho_vec){
  kmax <- length(rho_vec)
  kmax_use <- min(kmax, n-1)
  ks <- 1:kmax_use
  denom <- 1 + 2 * sum( (1 - ks/n) * abs(rho_vec[ks]) )
  n_eff <- n / denom
  return(n_eff)
}