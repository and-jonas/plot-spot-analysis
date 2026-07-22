
rm(list = ls())

.libPaths(c("~/R/library", .libPaths()))

# install required packages
list.of.packages <- c("tidyverse", "nlme", "pbmcapply")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, lib = "~/R/library", dependencies = TRUE, repos='https://stat.ethz.ch/CRAN/')

library(nlme)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(forcats)
library(tibble)
library(pbmcapply)

setwd("/agroscope/Data-Work-CH/22_Plant_Production-CH/224_Digitalisation/Jonas_Anderegg_Files/C_Manuscripts/2025_FocalStack_Downstream")

source("R/utils.R")
sub <- read.csv("data/subset.csv")

cat("Available cores:", parallel::detectCores(), "\n")
cat("Allocated cores:", Sys.getenv("SLURM_CPUS_PER_TASK"), "\n")

Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1
)

# ============================================================================== -
# 6) model-based; STEP 1 ----
# Auto-correlation in focal stacks
# ============================================================================== -

# >> 1) PLACL ----

# prepare data for modelling:
# to long format
mdat_stack <- sub %>% 
  pivot_wider(names_from = trait, values_from = value)
# logit-transformation
mdat_stack$placl_logit <- log((mdat_stack$PLACL + 0.001)/(1 - mdat_stack$PLACL + 0.001))
# convert to factors
mdat_stack <- mdat_stack %>%
  mutate(
    plot_id = factor(plot_id),
    position = factor(position),
    leaf_layer = factor(leaf_layer)
  )
# Add a unique identifier for each stack
mdat_stack$stack_uid = interaction(mdat_stack$plot_id, mdat_stack$position, mdat_stack$leaf_layer)

# bootstrap
nboot = 3
nsim = 3
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

bootstrap_results <- pbmclapply(
  1:nboot,
  function(b){
    
    cat("Bootstrap",b,"\n")
    
    # ---------------------------------------------------------- -
    # sample plots with replacement
    # ---------------------------------------------------------- -
    
    plots <- unique(mdat_stack$plot_id)
    sampled_plots <- sample(
      plots,
      length(plots),
      replace=TRUE
    )
    
    # make sure each sampled plot keeps a unique plot_id, even if the same plots
    # are sampled multiple times
    boot_data <- map2_dfr(
      sampled_plots,
      seq_along(sampled_plots),
      function(p, boot_id){
        mdat_stack %>%
          filter(plot_id == p) %>%
          mutate(
            plot_id = paste0(p, "_boot", boot_id)
          )
      }
    )
    # adjust the stack_uid
    boot_data <- boot_data %>%
      mutate(
        stack_uid = interaction(
          plot_id,
          position,
          leaf_layer
        )
      )
    
    # need unique factor levels
    boot_data <- droplevels(boot_data)
    
    # ---------------------------------------------------------- -
    # fit model
    # ---------------------------------------------------------- -
    
    # Fit the LMM
    model_boot <- try(
      lme(
        placl_logit ~ plot_id * leaf_layer,
        random =~1 | plot_id/position/leaf_layer,
        correlation = corAR1(form = ~ stack_image_id |plot_id/position/leaf_layer),
        data=boot_data,
        method="REML",
        na.action=na.exclude
      ),
      silent=TRUE
    )
    if(inherits(model_boot, "try-error")){
      return(NULL)
    }
    
    # ---------------------------------------------------------- -
    # simulate original scale ACF
    # ---------------------------------------------------------- -
    
    # get autocorrelation coefficient
    # this coefficient is valid for the logit-transformed data
    # logit-transformation 'stretches' values near 0 and 1
    # so a difference between PLACL 0 and 0.01 is 'amplified' compared to a difference between 
    # PLACL of 0.49 and 0.5; however, as breeders we would are more about values on the original PLACL scale ...
    phi1_stack_placl <- coef(model_boot$modelStruct$corStruct, unconstrained = FALSE)
    
    # Simulate and back-transform to estimate autocorrelation on original scale
    # model pieces
    sigma_resid <- sqrt(model_boot$sigma^2)
    fixef_vals <- fitted(model_boot, level = 0)
    # get random effects by stack
    ranef_df <- ranef(model_boot)
    
    # data layout for simulation: list of stacks with their row indices
    boot_data <- boot_data %>% mutate(stack_uid = fct_drop(stack_uid))
    stack_groups <- split(1:nrow(boot_data), boot_data$stack_uid)
    
    acf_sims <- matrix(NA, nrow = nsim, ncol = 10)  # lags 0..9
    for(sim in 1:nsim){
      print(sim)
      sim_vals_orig <- numeric(nrow(boot_data))
      for(si in seq_along(stack_groups)){
        rows <- stack_groups[[si]]
        n <- length(rows)
        # Build model-predicted mean on transformed scale for these rows
        # extract random effects per level
        identifiers <- str_split(names(stack_groups[si]), "\\.") %>% unlist()
        plot_id <- identifiers[1]
        position <- paste(plot_id, identifiers[2], sep = "/")
        leaf_layer <- paste(position, identifiers[3], sep = "/")
        ranef_plot <- ranef_df$plot_id
        ranef_pos  <- ranef_df$position
        ranef_leaf <- ranef_df$leaf_layer
        # build mu_trans by summing the random intercepts at each level
        mu_trans <- fixef_vals[rows] +
          unname(ranef_plot[plot_id, "(Intercept)"]) +
          unname(ranef_pos[position, "(Intercept)"]) +
          unname(ranef_leaf[leaf_layer, "(Intercept)"])
        # correlation matrix for this stack
        Cor <- get_cor_matrix(model_boot$modelStruct$corStruct, n)
        # we multiply the correlation matrix by the residual variance to get the full covariance matrix
        # resulting sigma encodes both variance and autocorrelation for the residuals within this stack
        Sigma <- (model_boot$sigma^2) * Cor
        # simulate residuals
        eps <- mvrnorm(1, mu = rep(0,n), Sigma = Sigma)
        # add the simulated residuals to the predicted means
        sim_trans <- mu_trans + eps
        # back-transform to the original scale
        # i.e., revert the logit transformation applied to the raw data
        sim_orig <- inv_logit_adjusted(sim_trans, eps = 0.001)
        sim_vals_orig[rows] <- sim_orig
      }
      # compute mean ACF across stacks for this simulated data set
      # compute per-stack acf and average (skip lag 0)
      acf_per_stack <- sapply(stack_groups, function(rows){
        r <- sim_vals_orig[rows]
        if(length(unique(r))>1) acf(r, lag.max=9, plot=FALSE, na.action=na.omit)$acf else rep(NA,9)
      })
      acf_per_stack_mat <- do.call(rbind, acf_per_stack)  # bind vectors row-wise
      acf_sims[sim, ] <- colMeans(acf_per_stack_mat, na.rm = TRUE)
    }
    
    # summarize
    model_based_acf_mean <- colMeans(acf_sims, na.rm=TRUE)
    model_based_acf_ci <- apply(acf_sims, 2, quantile, probs=c(0.025,0.975), na.rm=TRUE)
    model_based_acf_stack_placl <- tibble(lag = c(0:9), 
                                          trait = "PLACL",
                                          level = "stack", 
                                          mean = model_based_acf_mean,
                                          upper = model_based_acf_ci[1, ],
                                          lower = model_based_acf_ci[2, ])
    
    model_based_acf_stack <- bind_rows(model_based_acf_stack_placl)
    model_based_acf <- bind_rows(model_based_acf_stack)
    mbased <- model_based_acf
    # get point estimate
    n_eff_point_mbased <- mbased %>% group_by(trait, level) %>% nest() %>% 
      mutate(n = purrr::map_dbl(data, nrow),
             n_eff_point = purrr::map_dbl(data, ~n_eff_from_acf(nrow(.), .$mean)),
             n_eff_point_norm = n_eff_point/n)
    n_eff = n_eff_point_mbased$n_eff_point
    return(
      list(
        n_eff = n_eff,
        phi = phi1_stack_placl
      )
    )},
  mc.cores = ncores
)

saveRDS(
  bootstrap_results,
  file = "data/bootstrap_results.rds"
)

cat("Bootstrap finished successfully\n")