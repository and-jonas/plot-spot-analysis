
rm(list = ls())

start_time <- Sys.time()

.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1
)

# install required packages
list.of.packages <- c("nlme", "dplyr", "tidyr", "purrr", "stringr", "forcats", "tibble", "parallel", "MASS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, lib = "~/R/library", dependencies = TRUE, repos='https://stat.ethz.ch/CRAN/')

library(nlme)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(forcats)
library(tibble)
library(parallel)
library(MASS)

setwd("/agroscope/Data-Work-CH/22_Plant_Production-CH/224_Digitalisation/Jonas_Anderegg_Files/C_Manuscripts/2025_FocalStack_Downstream")
# setwd("O:/Data-Work/22_Plant_Production-CH/224_Digitalisation/Jonas_Anderegg_Files/C_Manuscripts/2025_FocalStack_Downstream")

source("R/utils.R")
sub <- read.csv("data/subset.csv")

cat("Available cores:", parallel::detectCores(), "\n")

# ============================================================================== -
# > 1) PLACL ----
# ============================================================================== -

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

# ---------------------------------------------------------- -
# Bootstrap
# ---------------------------------------------------------- -

# global parameters
nboot = 1000
nsim = 20
cat("Allocated cores:", Sys.getenv("SLURM_CPUS_PER_TASK"), "\n")
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))-2
cat("Using", ncores, "\n")
# ncores <- 1
# ncores <- 4

# set up cluster
cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  .libPaths(c("~/R/library", .libPaths()))
})

clusterEvalQ(cl, {
  library(nlme)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(tibble)
  library(MASS)
})

clusterExport(
  cl,
  c(
    "mdat_stack",
    "get_cor_matrix",
    "inv_logit_adjusted",
    "n_eff_from_acf",
    "nsim"
  )
)

bootstrap_results <- parLapply(
  cl, 
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
    
    # ========================================================================== -
    # >> a) focal planes ----
    # ========================================================================== -
    
    # ---------------------------------------------------------- -
    # fit model
    # ---------------------------------------------------------- -
    
    model_boot <- tryCatch(
      {
        lme(
          placl_logit ~ plot_id * leaf_layer,
          random =~1 | plot_id/position/leaf_layer,
          correlation = corAR1(form = ~ stack_image_id |plot_id/position/leaf_layer),
          data=boot_data,
          method="REML",
          na.action=na.exclude
        )
      },
      error = function(e) {
        message(
          "Bootstrap ", b,
          " failed: ",
          e$message
        )
        return(NULL)
      }
    )
    
    if(is.null(model_boot)){
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
      max_lag <- max(lengths(acf_per_stack))
      acf_per_stack_mat <- do.call(
        rbind,
        lapply(acf_per_stack, function(x) {
          length(x) <- max_lag   # extends with NA
          x
        })
      )
      acf_sims[sim, ] <- colMeans(acf_per_stack_mat, na.rm = TRUE)
    }
    
    # ---------------------------------------------------------- -
    # summarize and extract n_eff
    # ---------------------------------------------------------- -¨
    
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
    n_eff_foc_placl = n_eff_point_mbased$n_eff_point
    
    # ========================================================================== -
    # >> b) plot positions ----
    # ========================================================================== -

    # we obtain an estimate from the model in STEP 1 for each stack
    # Get the row indices actually used by the model, subset data
    used_rows <- as.numeric(rownames(model_boot[["residuals"]]))
    mdat_pos <- boot_data[used_rows, ]

    # aggregate fitted values by stack
    mdat_pos <- aggregate(fitted(model_boot),
                          by = list(stack_id = mdat_pos$stack_uid,
                                    plot_id = mdat_pos$plot_id,
                                    position = mdat_pos$position,
                                    leaf_layer = mdat_pos$leaf_layer),
                          FUN = mean)

    mdat_pos <- mdat_pos %>%
      dplyr::select(plot_id, position, leaf_layer, stack_id, x) %>%
      as_tibble() %>%
      arrange(plot_id, position, leaf_layer)
    # add a unique plot_layer identifier
    # here, specifying each plot-layer combo as random effect seems appropriate
    # because plot and leaf_layer can be considered crossed factors, rather than nested
    mdat_pos$plot_layer <- interaction(mdat_pos$plot_id, mdat_pos$leaf_layer, drop = TRUE)

    # Fit AR1 along positions within plots
    model_position_placl <- tryCatch(
      {
        lme(
          x ~ plot_id * leaf_layer,
          random = ~1 | plot_id/leaf_layer,
          correlation = corAR1(form = ~ as.numeric(position) | plot_id/leaf_layer),
          data = mdat_pos
        )
      },
      error = function(e) {
        message(
          "Bootstrap ", b,
          " failed: ",
          e$message
        )
        return(NULL)
      }
    )
    if(is.null(model_position_placl)){
      return(NULL)
    }
    
    # get autocorrelation coefficient
    # Valid on logit-transformed scale
    phi1_pos_placl <- coef(model_position_placl$modelStruct$corStruct, unconstrained = FALSE)

    # model pieces
    sigma_resid <- sqrt(model_position_placl$sigma^2)  # residual sd
    fixef_vals <- fitted(model_position_placl, level = 0)  # fixed-effect fitted values on transformed scale
    # get random effects by stack if you want to simulate them too:
    ranef_df <- ranef(model_position_placl)  # random intercepts per plot x layer

    # data layout for simulation: list of stacks with their row indices
    mdat_pos <- mdat_pos %>% mutate(plot_layer = fct_drop(plot_layer))
    stack_groups <- split(1:nrow(mdat_pos), mdat_pos$plot_layer)

    acf_sims <- matrix(NA, nrow = nsim, ncol = 25)  # lags 0..9
    for(sim in 1:nsim){
      print(sim)
      sim_vals_orig <- numeric(nrow(mdat_pos))
      for(si in seq_along(stack_groups)){
        rows <- stack_groups[[si]]
        n <- length(rows)
        # Build model-predicted mean on transformed scale for these rows
        # extract random effects per level
        identifiers <- str_split(names(stack_groups[si]), "\\.") %>% unlist()
        plot_id <- identifiers[1]
        leaf_layer <- paste(plot_id, identifiers[2], sep = "/")
        ranef_plot <- ranef_df$plot_id
        ranef_leaf <- ranef_df$leaf_layer
        # build mu_trans by summing the random intercepts at each level
        mu_trans <- fixef_vals[rows] +
          unname(ranef_plot[plot_id, "(Intercept)"]) +
          unname(ranef_leaf[leaf_layer, "(Intercept)"])
        # # Build model-predicted mean on transformed scale for these rows
        # mu_trans <- fixef_vals[rows] + unname(ranef_df[as.character(mdat_pos$plot_layer[rows]), "(Intercept)"])
        # correlation matrix for this stack
        Cor <- get_cor_matrix(model_position_placl$modelStruct$corStruct, n)
        # we multiply the correlation matrix by the residual variance to get the full covariance matrix
        # resulting sigma encodes both variance and autocorrelation for the residuals within this stack
        Sigma <- (model_position_placl$sigma^2) * Cor
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
        if(length(unique(r))>1) acf(r, lag.max=24, plot=FALSE, na.action=na.omit)$acf else rep(NA,25)
      })
      max_lag <- max(lengths(acf_per_stack))
      acf_per_stack_mat <- do.call(
        rbind,
        lapply(acf_per_stack, function(x) {
          length(x) <- max_lag   # extends with NA
          x
        })
      )
      acf_sims[sim, ] <- colMeans(acf_per_stack_mat, na.rm = TRUE)
    }

    # summarize
    model_based_acf_mean <- colMeans(acf_sims, na.rm=TRUE)
    model_based_acf_ci <- apply(acf_sims, 2, quantile, probs=c(0.025,0.975), na.rm=TRUE)
    model_based_acf_pos_placl <- tibble(lag = c(0:24),
                                        trait = "PLACL",
                                        level = "pos",
                                        mean = model_based_acf_mean,
                                        upper = model_based_acf_ci[1, ],
                                        lower = model_based_acf_ci[2, ])

    # model_based_acf_pos <- bind_rows(model_based_acf_pos_placl, model_based_acf_pos_rust)
    model_based_acf_pos <- bind_rows(model_based_acf_pos_placl)
    # model_based_acf <- bind_rows(model_based_acf_stack, model_based_acf_pos)
    model_based_acf <- bind_rows(model_based_acf_pos)
    mbased <- model_based_acf

    # get point estimate
    n_eff_point_mbased <- mbased %>% group_by(trait, level) %>% nest() %>%
      mutate(n = purrr::map_dbl(data, nrow),
             n_eff_point = purrr::map_dbl(data, ~n_eff_from_acf(nrow(.), .$mean)),
             n_eff_point_norm = n_eff_point/n)
    n_eff_pos_placl = n_eff_point_mbased$n_eff_point

    return(
      list(
        n_eff_foc_placl = n_eff_foc_placl,
        n_eff_pos_placl = n_eff_pos_placl,
        phi1_stack_placl = phi1_stack_placl,
        phi1_pos_placl = phi1_pos_placl
      )
    )
    }
)

stopCluster(cl)

saveRDS(
  bootstrap_results,
  file = "data/bootstrap_results_PLACL.rds"
)

end_time <- Sys.time()

cat(
  "Total runtime:",
  round(difftime(end_time, start_time, units = "mins"), 2),
  "minutes\n"
)

# ============================================================================== -
# > 2) Rust ----
# ============================================================================== -

# prepare data for modelling:
# to long format
mdat_stack <- sub %>%
  pivot_wider(names_from = trait, values_from = value)
# add a small constant to avoid -inf in log(x)
mdat_stack$Pustules <- mdat_stack$Pustules + 0.1
# log-transformation
mdat_stack$rust_log <- log(mdat_stack$Pustules)
# convert to factors
mdat_stack <- mdat_stack %>%
  mutate(
    plot_id = factor(plot_id),
    position = factor(position),
    leaf_layer = factor(leaf_layer)
  )
# Add a unique identifier for each stack
mdat_stack$stack_uid = interaction(mdat_stack$plot_id, mdat_stack$position, mdat_stack$leaf_layer)

# ---------------------------------------------------------- -
# Bootstrap
# ---------------------------------------------------------- -

# parameters as for PLACL
#
# set up cluster
cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  .libPaths(c("~/R/library", .libPaths()))
})
clusterEvalQ(cl, {
  library(nlme)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(tibble)
  library(MASS)
})
clusterExport(
  cl,
  c(
    "mdat_stack",
    "get_cor_matrix",
    "inv_logit_adjusted",
    "n_eff_from_acf",
    "nsim"
  )
)

bootstrap_results <- parLapply(
  cl,
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

    # ========================================================================== -
    # >> a) focal planes ----
    # ========================================================================== -

    # ---------------------------------------------------------- -
    # fit model
    # ---------------------------------------------------------- -

    # Fit the LMM
    model_boot <- tryCatch(
      {
        lme(
          rust_log ~ plot_id * leaf_layer,
          random = ~ 1 | plot_id/position/leaf_layer,
          correlation = corAR1(form = ~ stack_image_id | plot_id/position/leaf_layer),
          data = boot_data,
          method = "REML",
          na.action=na.exclude
        )
      },
      error = function(e) {
        message(
          "Bootstrap ", b,
          " failed: ",
          e$message
        )
        return(NULL)
      }
    )

    if(is.null(model_boot)){
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
    phi1_stack_rust <- coef(model_boot$modelStruct$corStruct, unconstrained = FALSE)

    # SIMULATION
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
        sim_orig <- inv_log_adjusted(sim_trans, eps = 0.1)
        sim_vals_orig[rows] <- sim_orig
      }
      # compute mean ACF across stacks for this simulated data set
      # compute per-stack acf and average (skip lag 0)
      acf_per_stack <- sapply(stack_groups, function(rows){
        r <- sim_vals_orig[rows]
        if(length(unique(r))>1) acf(r, lag.max=9, plot=FALSE, na.action=na.omit)$acf else rep(NA,9)
      })
      max_lag <- max(lengths(acf_per_stack))
      acf_per_stack_mat <- do.call(
        rbind,
        lapply(acf_per_stack, function(x) {
          length(x) <- max_lag   # extends with NA
          x
        })
      )      acf_sims[sim, ] <- colMeans(acf_per_stack_mat, na.rm = TRUE)
    }

    # summarize
    model_based_acf_mean <- colMeans(acf_sims, na.rm=TRUE)
    model_based_acf_ci <- apply(acf_sims, 2, quantile, probs=c(0.025,0.975), na.rm=TRUE)
    model_based_acf_stack_rust <- tibble(lag = c(0:9),
                                         trait = "Pustules",
                                         level = "stack",
                                         mean = model_based_acf_mean,
                                         upper = model_based_acf_ci[1, ],
                                         lower = model_based_acf_ci[2, ])

    model_based_acf_stack <- bind_rows(model_based_acf_stack_rust)
    model_based_acf <- bind_rows(model_based_acf_stack)
    mbased <- model_based_acf
    # get point estimate
    n_eff_point_mbased <- mbased %>% group_by(trait, level) %>% nest() %>%
      mutate(n = purrr::map_dbl(data, nrow),
             n_eff_point = purrr::map_dbl(data, ~n_eff_from_acf(nrow(.), .$mean)),
             n_eff_point_norm = n_eff_point/n)
    n_eff_foc_rust <- n_eff_point_mbased$n_eff_point

    # ========================================================================== -
    # >> b) plot positions ----
    # ========================================================================== -

    # we obtain an estimate from the model in STEP 1 for each stack
    # Get the row indices actually used by the model, subset data
    used_rows <- as.numeric(rownames(model_boot[["residuals"]]))
    mdat_pos <- boot_data[used_rows, ]

    # aggregate fitted values by stack
    mdat_pos <- aggregate(fitted(model_boot),
                          by = list(stack_id = mdat_pos$stack_uid,
                                    plot_id = mdat_pos$plot_id,
                                    position = mdat_pos$position,
                                    leaf_layer = mdat_pos$leaf_layer),
                          FUN = mean)

    mdat_pos <- mdat_pos %>%
      dplyr::select(plot_id, position, leaf_layer, stack_id, x) %>%
      as_tibble() %>%
      arrange(plot_id, position, leaf_layer)
    # add a unique plot_layer identifier
    # here, specifying each plot-layer combo as random effect seems appropriate
    # because plot and leaf_layer can be considered crossed factors, rather than nested
    mdat_pos$plot_layer <- interaction(mdat_pos$plot_id, mdat_pos$leaf_layer, drop = TRUE)

    # Fit AR1 along positions within plots
    model_position_rust <- lme(
      x ~ plot_id * leaf_layer,
      random = ~1 | plot_id/leaf_layer,
      correlation = corAR1(form = ~ as.numeric(position) | plot_id/leaf_layer),
      data = mdat_pos
    )

    # get autocorrelation coefficient
    # Valid on logit-transformed scale
    phi1_pos_rust <- coef(model_position_rust$modelStruct$corStruct, unconstrained = FALSE)

    # model pieces
    sigma_resid <- sqrt(model_position_rust$sigma^2)  # residual sd
    fixef_vals <- fitted(model_position_rust, level = 0)  # fixed-effect fitted values on transformed scale
    # get random effects by stack if you want to simulate them too:
    ranef_df <- ranef(model_position_rust)  # random intercepts per plot x layer

    # data layout for simulation: list of stacks with their row indices
    mdat_pos <- mdat_pos %>% mutate(plot_layer = fct_drop(plot_layer))
    stack_groups <- split(1:nrow(mdat_pos), mdat_pos$plot_layer)

    acf_sims <- matrix(NA, nrow = nsim, ncol = 25)  # lags 0..9
    for(sim in 1:nsim){
      print(sim)
      sim_vals_orig <- numeric(nrow(mdat_pos))
      for(si in seq_along(stack_groups)){
        rows <- stack_groups[[si]]
        n <- length(rows)
        # Build model-predicted mean on transformed scale for these rows
        # extract random effects per level
        identifiers <- str_split(names(stack_groups[si]), "\\.") %>% unlist()
        plot_id <- identifiers[1]
        leaf_layer <- paste(plot_id, identifiers[2], sep = "/")
        ranef_plot <- ranef_df$plot_id
        ranef_leaf <- ranef_df$leaf_layer
        # build mu_trans by summing the random intercepts at each level
        mu_trans <- fixef_vals[rows] +
          unname(ranef_plot[plot_id, "(Intercept)"]) +
          unname(ranef_leaf[leaf_layer, "(Intercept)"])
        # # Build model-predicted mean on transformed scale for these rows
        # mu_trans <- fixef_vals[rows] + unname(ranef_df[as.character(mdat_pos$plot_layer[rows]), "(Intercept)"])
        # correlation matrix for this stack
        Cor <- get_cor_matrix(model_position_rust$modelStruct$corStruct, n)
        # we multiply the correlation matrix by the residual variance to get the full covariance matrix
        # resulting sigma encodes both variance and autocorrelation for the residuals within this stack
        Sigma <- (model_position_rust$sigma^2) * Cor
        # simulate residuals
        eps <- mvrnorm(1, mu = rep(0,n), Sigma = Sigma)
        # add the simulated residuals to the predicted means
        sim_trans <- mu_trans + eps
        # back-transform to the original scale
        # i.e., revert the logit transformation applied to the raw data
        sim_orig <- inv_log_adjusted(sim_trans, eps = 0.001)
        sim_vals_orig[rows] <- sim_orig
      }
      # compute mean ACF across stacks for this simulated data set
      # compute per-stack acf and average (skip lag 0)
      acf_per_stack <- sapply(stack_groups, function(rows){
        r <- sim_vals_orig[rows]
        if(length(unique(r))>1) acf(r, lag.max=24, plot=FALSE, na.action=na.omit)$acf else rep(NA,25)
      })
      max_lag <- max(lengths(acf_per_stack))
      acf_per_stack_mat <- do.call(
        rbind,
        lapply(acf_per_stack, function(x) {
          length(x) <- max_lag   # extends with NA
          x
        })
      )
      acf_sims[sim, ] <- colMeans(acf_per_stack_mat, na.rm = TRUE)
    }

    # summarize
    model_based_acf_mean <- colMeans(acf_sims, na.rm=TRUE)
    model_based_acf_ci <- apply(acf_sims, 2, quantile, probs=c(0.025,0.975), na.rm=TRUE)
    model_based_acf_pos_rust <- tibble(lag = c(0:24),
                                       trait = "Pustules",
                                       level = "pos",
                                       mean = model_based_acf_mean,
                                       upper = model_based_acf_ci[1, ],
                                       lower = model_based_acf_ci[2, ])

    # model_based_acf_pos <- bind_rows(model_based_acf_pos_placl, model_based_acf_pos_rust)
    model_based_acf_pos <- bind_rows(model_based_acf_pos_rust)
    # model_based_acf <- bind_rows(model_based_acf_stack, model_based_acf_pos)
    model_based_acf <- bind_rows(model_based_acf_pos)
    mbased <- model_based_acf

    # get point estimate
    n_eff_point_mbased <- mbased %>% group_by(trait, level) %>% nest() %>%
      mutate(n = purrr::map_dbl(data, nrow),
             n_eff_point = purrr::map_dbl(data, ~n_eff_from_acf(nrow(.), .$mean)),
             n_eff_point_norm = n_eff_point/n)
    n_eff_pos_rust = n_eff_point_mbased$n_eff_point

    return(
      list(
        n_eff_foc_rust = n_eff_foc_rust,
        n_eff_pos_rust = n_eff_pos_rust,
        phi1_stack_rust = phi1_stack_rust,
        phi1_pos_rust = phi1_pos_rust
      )
    )}
)

stopCluster(cl)

saveRDS(
  bootstrap_results,
  file = "data/bootstrap_results_rust.rds"
)

end_time <- Sys.time()

cat(
  "Total runtime:",
  round(difftime(end_time, start_time, units = "mins"), 2),
  "minutes\n"
)
