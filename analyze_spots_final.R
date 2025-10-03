
# ============================================================================== -
# HEADER ----
# ============================================================================== -

rm(list = ls())

library(nlme)
library(MASS)
library(tidyverse)

setwd("O:/Data-Work/22_Plant_Production-CH/224_Digitalisation/Jonas_Anderegg_Files/C_Manuscripts/2025_FocalStack")

# basic plot theme
base_theme <- theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# helpers
# invert logit
inv_logit_adjusted <- function(x, eps = 0.001) {
  # Standard inverse logit
  p <- exp(x) / (1 + exp(x))
  # Undo the + eps adjustment approximately
  (p - eps) / (1 - 2*eps)
}
# invert log
inv_log_adjusted <- function(x, eps = 0.1) {
  # Undo log(x + eps) -> return x
  orig <- exp(x) - eps
  # Clamp to >= 0 to avoid negative values
  orig[orig < 0] <- 0
  return(orig)
}
# get full correlation matrix
get_cor_matrix <- function(corStruct, n){
  phi <- coef(corStruct, unconstrained = FALSE)
  mat <- phi ^ (abs(outer(1:n,1:n, "-")))
  return(mat)
}

# ============================================================================== -
# 0) pre-process data ----
# Minor issues: 
# - few multiple stacks per position x layer x plot
# - few incomplete stacks
# ============================================================================== -

# read raw data
data_raw <- read_csv("overall_canopy_results.csv")

# check if there is one stack per group, and if there are 10 images per stack

check_unique <- data_raw %>% 
  # subset biggest campaign
  dplyr::filter(date %in% c("2023-06-06", "2023-06-07")) %>%
  # add a camera id
  mutate(cam_id = str_sub(imagefile_uid, 1, 3)) %>%
  # get number of images per stack
  group_by(plot_id, leaf_layer, position, cam_id) %>%
  nest() %>%
  mutate(n = purrr::map_dbl(data, nrow))

# identify problems
not_complete_not_unique <- check_unique %>%
  filter(n != 10)
# does not select a unique stack for all plot x position (some have two stacks)
# some stacks for plot 209 are incomplete (9 and 1 images)

# keep only complete stacks (dropping 2 incomplete stacks)
complete_not_unique <- check_unique %>% 
  filter(n %in% c(10, 20))

# split stacks where multiple exist per plot
data <- complete_not_unique %>% 
  # remove incomplete stacks
  mutate(data = map(data, ~ {
    n <- nrow(.x)
    if (n == 10) {
      list(.x)  # keep as single list element
    } else if (n == 20) {
      print("splitting")
      split(.x, rep(1:2, each = 10))  # split into two tibbles
    } else {
      stop("Unexpected stack size: ", n)  # sanity check
    }
  })) %>%
  unnest_longer(data) %>%
  group_by(plot_id, leaf_layer, position) %>%
  mutate(stack_id = row_number()) %>%
  ungroup() %>% 
  dplyr::select(-data_id, -n) %>% unnest(cols = c(data)) %>% ungroup()

# check success
data %>% 
  group_by(plot_id, leaf_layer, position, stack_id) %>%
  nest() %>%
  mutate(n = purrr::map_dbl(data, nrow)) %>% pull(n) %>% unique()
# OK

# still, multiple stacks per plot x layer x position cause complications 
# with plotting and acf calculations
# --> I just drop those plot x layer x position combos
data <- data %>% 
  filter(stack_id == 1)

# ============================================================================== -
# Add design
# ============================================================================== -

# add leaf angle to design info
ang <- c(rep("L", 6), rep("H", 6))
gen <- c("GULLIVER", "KOBRA PLUS", "HISTORY", "RAISON", "BRANDO", "POTENZIAL",
         "KWS SCIROCCO LP 509.3.04", "XENOS", "TAIFUN", "DIAVEL", "RETRO", "BORNEO")
info <- tibble(ang, gen)

# get experimental design
design <- data_raw %>% dplyr::filter(date %in% c("2023-06-06", "2023-06-07")) %>% 
  full_join(., info, by = c("genotype_name" = "gen")) %>% 
  dplyr::select(plot_id, ang, genotype_name) %>% unique()

# select data corresponding to the largest campaign
# rearrange data
sub <- data %>% 
  dplyr::select(placl, `rust_density_1e-6`, plot_id, leaf_layer, position, stack_image_id) %>% 
  pivot_longer(cols = c(placl, `rust_density_1e-6`), names_to = "trait") %>% 
  # add readable trait names
  mutate(trait = ifelse(trait == "placl", "PLACL", "Pustules")) %>% 
  #add experimental design
  full_join(., design) %>% 
  # remove rows with missing trait values (NA because nothing was in focus)
  filter(!is.na(value)) %>% droplevels(.)

# trait values per genotype
ggplot(sub) +
  geom_boxplot(aes(x = genotype_name, y = value)) + 
  facet_wrap(~trait, scales = "free")

# ============================================================================== -
# 1) prelim: stackpos vs. placl ----
# ============================================================================== -

# visualize position-in-stack vs. trait value relationship
# position-in-stack is associated with trait value!
# positions in the beginning of the stack tend to have lower trait values 
# than positions further down the stack
pos_eff <- ggplot(sub) +
  geom_density(aes(x = value, color = as.factor(stack_image_id)), alpha = 0.4) +
  theme(panel.grid.minor = element_blank())+
  facet_wrap(~trait, scales = "free_y")
png("stackpos_traits_all.png", width = 6, height = 6, units = 'in', res = 400)
plot(pos_eff)
dev.off()

# plot for each layer
# in absolute terms, this is particularly pronounced for leaf layer 3, less for 2, and even less for 1
pdat <- sub %>% 
  pivot_wider(values_from = value, names_from = trait) %>% 
  # filter out 'extreme' observations to improve visibility of trends in density plots
  filter(!(leaf_layer == 1 & PLACL > 0.03)) %>% 
  filter(!(leaf_layer == 2 & PLACL > 0.25)) %>%
  filter(!(Pustules > 5)) %>% 
  pivot_longer(cols = c(PLACL, Pustules), names_to = "trait", values_to = "value")

pdat <- pdat %>%
  mutate(facet_label = paste0("Layer ", leaf_layer, " - ", trait))

stack_density <- ggplot(pdat) + 
  stat_density(
    aes(x = value, color = as.factor(stack_image_id)),
    geom = "line", 
    position = "identity",
    size = 0.9
  ) +
  facet_wrap(
    ~facet_label, 
    ncol = 2, 
    scales = "free"
  ) +
  scale_color_viridis_d(name = "Image in stack", option = "plasma") +
  guides(
    color = guide_legend(
      override.aes = list(size = 1),
      nrow = 4
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = c(0.25, 0.9),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7, face = "bold")
  ) +
  labs(
    x = "Trait value",   
    y = "Density"
  )

png("stackpos_traits_layers.png", width = 5, height = 7, units = 'in', res = 400)
plot(stack_density)
dev.off()

# ============================================================================== -
# 2) acf: focal stacks ----
# ============================================================================== -

# estimate acf on cleaned data
acf_list_clean <- sub %>%
  group_by(plot_id, leaf_layer, position, trait) %>% 
  summarise(acf_vals = list(acf(value, plot = FALSE, lag.max = 9, na.action = na.pass)$acf),
            .groups = "drop") %>% 
  full_join(., design)

# convert to long format tibble for plotting
acf_long <- acf_list_clean %>%
  mutate(acf_vals = map(acf_vals, as.numeric)) %>%
  unnest_longer(acf_vals) %>%
  group_by(plot_id, leaf_layer, position, trait) %>%
  mutate(lag = row_number() - 1) %>%
  ungroup()

# plot autocorrelation pattern for each individual stack
stacks <- ggplot(acf_long, aes(x = lag, y = acf_vals, group = interaction(plot_id, leaf_layer, position, trait))) +
  geom_line(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  labs(x = "Lag", y = "Autocorrelation", color = "Genotype", fill = "Genotype") +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  facet_wrap(~trait) +
  theme(panel.grid.minor = element_blank()) +
  base_theme
png("per_stack_traits.png", width = 8, height = 4, units = 'in', res = 400)
plot(stacks)
dev.off()

# ============================================================================== -
# 3) acf: focal stacks: mean ± SD ----
# ============================================================================== -

# get overall means and standard deviation
pdat_summary_stack <- acf_long %>%
  group_by(lag, trait) %>%
  summarise(
    mean_acf = mean(acf_vals, na.rm = TRUE),
    sd_acf = sd(acf_vals, na.rm = TRUE),
    .groups = "drop"
  )

# plot mean ± SD
stack_all <- ggplot(pdat_summary_stack, aes(x = lag, y = mean_acf, group = trait, color = trait)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = trait),
              alpha = 0.2, color = NA) +
  labs(x = "Lag", y = "Autocorrelation") +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  base_theme
  
png("acf_stack_all.png", width = 6, height = 6, units = 'in', res = 400)
plot(stack_all)
dev.off()

# ============================================================================== -
# 4) acf: positions ----
# ============================================================================== -

# get simple arithmetic stack means for this preliminary analysis
stack_mean <- sub %>% 
  group_by(plot_id, leaf_layer, position, trait) %>% 
  summarize(mean = mean(value, na.rm = T),
            .groups = "drop")

# rearrange data and estimate acf
acf_list_clean <- stack_mean %>%
  group_by(plot_id, leaf_layer, trait) %>% 
  nest() %>% 
  mutate(n_pos = purrr::map_dbl(data, nrow)) %>% 
  # number of positions varies from 18 to 72 (why?)
  # total of 45 plot-by-layer combos available (why not 16*3 = 48?)
  # filter out plots with number of positions too different from 32
  # physical distance between positions would be too different in those plots
  # this leaves 39 plot-layer combos for analysis
  filter(n_pos >=30 && n_pos <= 40) %>% 
  unnest(cols = c(data)) %>% 
  summarise(acf_vals = list(acf(mean, plot = FALSE, lag.max = 24, na.action = na.pass)$acf),
            .groups = "drop") %>% 
  full_join(., design)

# convert to long format tibble for plotting
acf_long <- acf_list_clean %>%
  mutate(acf_vals = map(acf_vals, as.numeric)) %>%
  unnest_longer(acf_vals) %>%
  group_by(plot_id, leaf_layer, trait) %>%
  mutate(lag = row_number() - 1) %>%
  ungroup()

# plot autocorrelation pattern for each individual stack
positions <- ggplot(acf_long, aes(x = lag, y = acf_vals, group = interaction(plot_id, leaf_layer, trait))) +
  geom_line(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  labs(x = "Lag", y = "Autocorrelation", color = "Genotype", fill = "Genotype") +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  facet_wrap(~trait) +
  base_theme
png("acf_pos_plot_layer.png", width = 8, height = 4, units = 'in', res = 400)
plot(positions)
dev.off()

# ============================================================================== -
# 5) acf: positions: mean ± SD ---- ----
# ============================================================================== -

# get overall means and standard deviation
pdat_summary_pos <- acf_long %>%
  group_by(lag, trait) %>%
  summarise(
    mean_acf = mean(acf_vals, na.rm = TRUE),
    sd_acf = sd(acf_vals, na.rm = TRUE),
    .groups = "drop"
  )

# plot mean ± SD
pos_all <- ggplot(pdat_summary_pos, aes(x = lag, y = mean_acf, group = trait, color = trait)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = trait),
              alpha = 0.2, color = NA) +
  labs(x = "Lag", y = "Autocorrelation") +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  base_theme
png("acf_pos_all.png", width = 3, height = 3, units = 'in', res = 400)
plot(pos_all)
dev.off()

# COMBINED PLOT (STACK AND POS)
pdat_summary_pos$level <- "pos"
pdat_summary_stack$level  <- "stack"
pdat_summary_both <- bind_rows(pdat_summary_pos, pdat_summary_stack)
both <- ggplot(pdat_summary_both, aes(x = lag, y = mean_acf, group = trait, color = trait)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0, color = "grey70") +
  geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = trait),
              alpha = 0.2, color = NA) +
  labs(x = "Lag", y = "Autocorrelation") +
  scale_color_manual(values = c("#5CA37C", "#A35C83")) +
  scale_fill_manual(values = c("#5CA37C", "#A35C83")) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  facet_wrap(~level, ncol = 1, scales = "free_x",
             labeller = labeller(level = c(
               "pos" = "Position in Plot",
               "stack" = "Position in Stack"
             ))) +
  base_theme

png("acf_both_all.png", width = 3, height = 3, units = 'in', res = 400)
plot(both)
dev.off()

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

# Fit an LMM
# we have "pseudo-replicates" from the positions (1...32) for each plot.
# Actually, we would need to specify 'position nested in plot' as a random factor to account for this. 
# However, it is apparently not possible to have different formulae for the random and the correlation term; 
# and nested terms are not allowed in the correlation term. Therefore, we include a random effect for each stack instead. 
# According to my understanding, I cannot have nested random effects for plot/position and an AR(1) correlation along images in stacks at the same time. 
# So this is not perfect, but OK for descriptive purposes? We don't do formal inference on plot-level effects, the main goal 
# is to quantify auto-correlation based on the full data set. 
# Leaf layer can be considered crossed with plot (leaf layer 1 has the same/a very similar meaning across all plots)
# We include leaf layer as an additional fixed effect

# stack_image_id is 1..10, thus identifying the "depth" of the image in each stack
# stack_uid is a unique identifier for each focal stack

# without plot x layer interaction
model_stacks_0 <- lme(
  placl_logit ~ plot_id + leaf_layer,
  random = ~1 | stack_uid,
  correlation = corAR1(form = ~ stack_image_id | stack_uid),
  data = mdat_stack,
  method = "REML",
  na.action=na.exclude
)
# with plot x layer interaction
model_stacks_1 <- lme(
  placl_logit ~ plot_id * leaf_layer,
  random = ~1 | stack_uid,
  correlation = corAR1(form = ~ stack_image_id | stack_uid),
  data = mdat_stack,
  method = "REML",
  na.action=na.exclude
)
AIC(model_stacks_0, model_stacks_1)  # preference for model 1
model_stacks_placl <- model_stacks_1 

plot(model_stacks_placl)
qqnorm(resid(model_stacks_placl))
qqline(resid(model_stacks_placl))
hist(resid(model_stacks_placl))
# OK, I guess? Best we can do?

# get autocorrelation coefficient
# this coefficient is valid for the logit-transformed data
# logit-transformation 'stretches' values near 0 and 1
# so a difference between PLACL 0 and 0.1 is 'amplified' compared to a difference between 
# PLACL of 0.49 and 0.5; however, as breeders we would are more about values on the original PLACL scale ...
phi1_stack_placl <- coef(model_stacks_placl$modelStruct$corStruct, unconstrained = FALSE)

# SIMULATION
# model pieces
sigma_resid <- sqrt(model_stacks_placl$sigma^2)  # residual sd
fixef_vals <- fitted(model_stacks_placl, level = 0)  # fixed-effect fitted values on transformed scale
# get random effects by stack if you want to simulate them too:
ranef_df <- ranef(model_stacks_placl)  # random intercepts per stack

# data layout for simulation: list of stacks with their row indices
mdat_stack <- mdat_stack %>% mutate(stack_uid = fct_drop(stack_uid))
stack_groups <- split(1:nrow(mdat_stack), mdat_stack$stack_uid)

nsim <- 20
acf_sims <- matrix(NA, nrow = nsim, ncol = 10)  # lags 0..9

for(sim in 1:nsim){
  print(sim)
  sim_vals_orig <- numeric(nrow(mdat_stack))
  for(si in seq_along(stack_groups)){
    rows <- stack_groups[[si]]
    n <- length(rows)
    # Build model-predicted mean on transformed scale for these rows
    mu_trans <- fixef_vals[rows] + unname(ranef_df[as.character(mdat_stack$stack_uid[rows]), "(Intercept)"])
    # correlation matrix for this stack
    Cor <- get_cor_matrix(model_stacks_placl$modelStruct$corStruct, n)
    # we multiply the correlation matrix by the residual variance to get the full covariance matrix
    # resulting sigma encodes both variance and autocorrelation for the residuals within this stack
    Sigma <- (model_stacks_placl$sigma^2) * Cor
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

# >> 2) Rust ----

# prepare data for modelling:
# to long format
mdat_stack <- sub %>% 
  pivot_wider(names_from = trait, values_from = value)
ggplot(mdat_stack) +
  geom_histogram(aes(x = Pustules), bins = 100) +
  scale_x_continuous(limits = c(-0.1, 10))
# add a small constant to avoid -inf in log(x)
mdat_stack$Pustules <- mdat_stack$Pustules + 0.1
ggplot(mdat_stack) +
  geom_histogram(aes(x = Pustules), bins = 100) +
  scale_x_continuous(limits = c(-0.1, 10))

# logit-transformation
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

# Fit an LMM
# we have "pseudo-replicates" from the positions (1...32) for each plot.
# Actually, we would need to specify 'position nested in plot' as a random factor to account for this. 
# However, it is apparently not possible to have different formulae for the random and the correlation term; 
# and nested terms are not allowed in the correlation term. Therefore, we include a random effect for each stack instead. 
# According to my understanding, I cannot have nested random effects for plot/position and an AR(1) correlation along images in stacks at the same time. 
# So this is not perfect, but OK for descriptive purposes? We don't do formal inference on plot-level effects, the main goal 
# is to quantify auto-correlation based on the full data set. 
# Leaf layer can be considered crossed with plot (leaf layer 1 has the same/a very similar meaning across all plots)
# We include leaf layer as an additional fixed effect

# stack_image_id is 1..10, thus identifying the "depth" of the image in each stack
# stack_uid is a unique identifier for each focal stack

# # without plot x layer interaction
# model_stacks_0 <- lme(
#   rust_log ~ plot_id + leaf_layer,
#   random = ~1 | stack_uid,
#   correlation = corAR1(form = ~ stack_image_id | stack_uid),
#   data = mdat_stack,
#   method = "REML",
#   na.action=na.exclude
# )
# with plot x layer interaction
model_stacks_rust <- lme(
  rust_log ~ plot_id * leaf_layer,
  random = ~1 | stack_uid,
  correlation = corAR1(form = ~ stack_image_id | stack_uid),
  data = mdat_stack,
  method = "REML",
  na.action=na.exclude
)
# AIC(model_stacks_0, model_stacks_1)  # preference for model 1
plot(model_stacks_rust)
qqnorm(resid(model_stacks_rust))
qqline(resid(model_stacks_rust))
hist(resid(model_stacks_rust))
# OK, I guess? Best we can do?

# get autocorrelation coefficient
# this coefficient is valid for the logit-transformed data
# logit-transformation 'stretches' values near 0 and 1
# so a difference between PLACL 0 and 0.1 is 'amplified' compared to a difference between 
# PLACL of 0.49 and 0.5; however, as breeders we would are more about values on the original PLACL scale ...
phi1_stack_rust <- coef(model_stacks_rust$modelStruct$corStruct, unconstrained = FALSE)

# SIMULATION
# model pieces
sigma_resid <- sqrt(model_stacks_rust$sigma^2)  # residual sd
fixef_vals <- fitted(model_stacks_rust, level = 0)  # fixed-effect fitted values on transformed scale
# get random effects by stack if you want to simulate them too:
ranef_df <- ranef(model_stacks_rust)  # random intercepts per stack

# data layout for simulation: list of stacks with their row indices
mdat_stack <- mdat_stack %>% mutate(stack_uid = fct_drop(stack_uid))
stack_groups <- split(1:nrow(mdat_stack), mdat_stack$stack_uid)

nsim <- 20
acf_sims <- matrix(NA, nrow = nsim, ncol = 10)  # lags 0..9

for(sim in 1:nsim){
  print(sim)
  sim_vals_orig <- numeric(nrow(mdat_stack))
  for(si in seq_along(stack_groups)){
    rows <- stack_groups[[si]]
    n <- length(rows)
    # Build model-predicted mean on transformed scale for these rows
    mu_trans <- fixef_vals[rows] + unname(ranef_df[as.character(mdat_stack$stack_uid[rows]), "(Intercept)"])
    # correlation matrix for this stack
    Cor <- get_cor_matrix(model_stacks_rust$modelStruct$corStruct, n)
    # we multiply the correlation matrix by the residual variance to get the full covariance matrix
    # resulting sigma encodes both variance and autocorrelation for the residuals within this stack
    Sigma <- (model_stacks_rust$sigma^2) * Cor
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
  acf_per_stack_mat <- do.call(rbind, acf_per_stack)  # bind vectors row-wise
  acf_sims[sim, ] <- colMeans(acf_per_stack_mat, na.rm = TRUE)
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

# bind results
model_based_acf_stack <- bind_rows(model_based_acf_stack_placl, model_based_acf_stack_rust)

# ============================================================================== -
# 7) model-based; STEP 2 ----
# Auto-correlation along positions-in-plot
# ============================================================================== -

# >> 1) PLACL ----

# we obtain an estimate from the model in STEP 1 for each stack
# Get the row indices actually used by the model, subset data
used_rows <- as.numeric(rownames(model_stacks_placl[["residuals"]]))
mdat_pos <- mdat_stack[used_rows, ]

# aggregate fitted values by stack
mdat_pos <- aggregate(fitted(model_stacks_placl), 
                        by = list(stack_id = mdat_pos$stack_uid,
                                   plot_id = mdat_pos$plot_id,
                                   position = mdat_pos$position,
                                   leaf_layer = mdat_pos$leaf_layer), 
                         FUN = mean)

# invert logit back to original scale
inv_logit <- function(x) exp(x) / (1 + exp(x))
mean_fitted_placl <- mean(inv_logit(mdat_pos$x))  # seems plausible

mdat_pos <- mdat_pos %>% 
  dplyr::select(plot_id, position, leaf_layer, stack_id, x) %>%
  as_tibble() %>% 
  arrange(plot_id, position, leaf_layer)
# add a unique plot_layer identifier
# here, specifying each plot-layer combo as random effect seems appropriate 
# because plot and leaf_layer can be considered crossed factors, rather than nested
mdat_pos$plot_layer <- interaction(mdat_pos$plot_id, mdat_pos$leaf_layer, drop = TRUE)

# Fit AR1 along positions within plots
model_position <- lme(
  x ~ plot_id * leaf_layer,
  random = ~1 | plot_layer,
  correlation = corAR1(form = ~ as.numeric(position) | plot_layer),
  data = mdat_pos
)
plot(model_position)
qqnorm(resid(model_position))
qqline(resid(model_position))
hist(resid(model_position))
# Seems OK.

# get autocorrelation coefficient
# Valid on logit-transformed scale
phi1_pos_placl <- coef(model_position$modelStruct$corStruct, unconstrained = FALSE)

# model pieces
sigma_resid <- sqrt(model_position$sigma^2)  # residual sd
fixef_vals <- fitted(model_position, level = 0)  # fixed-effect fitted values on transformed scale
# get random effects by stack if you want to simulate them too:
ranef_df <- ranef(model_position)  # random intercepts per plot x layer

# data layout for simulation: list of stacks with their row indices
mdat_pos <- mdat_pos %>% mutate(plot_layer = fct_drop(plot_layer))
stack_groups <- split(1:nrow(mdat_pos), mdat_pos$plot_layer)

nsim <- 20
acf_sims <- matrix(NA, nrow = nsim, ncol = 25)  # lags 0..9

for(sim in 1:nsim){
  print(sim)
  sim_vals_orig <- numeric(nrow(mdat_pos))
  for(si in seq_along(stack_groups)){
    rows <- stack_groups[[si]]
    n <- length(rows)
    # Build model-predicted mean on transformed scale for these rows
    mu_trans <- fixef_vals[rows] + unname(ranef_df[as.character(mdat_pos$plot_layer[rows]), "(Intercept)"])
    # correlation matrix for this stack
    Cor <- get_cor_matrix(model_position$modelStruct$corStruct, n)
    # we multiply the correlation matrix by the residual variance to get the full covariance matrix
    # resulting sigma encodes both variance and autocorrelation for the residuals within this stack
    Sigma <- (model_position$sigma^2) * Cor
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
  acf_per_stack_mat <- do.call(rbind, acf_per_stack)  # bind vectors row-wise
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

# >> 2) Rust ----

# we obtain an estimate from the model in STEP 1 for each stack
# Get the row indices actually used by the model, subset data
used_rows <- as.numeric(rownames(model_stacks_rust[["residuals"]]))
mdat_pos <- mdat_stack[used_rows, ]

# aggregate fitted values by stack
mdat_pos <- aggregate(fitted(model_stacks_rust), 
                      by = list(stack_id = mdat_pos$stack_uid,
                                plot_id = mdat_pos$plot_id,
                                position = mdat_pos$position,
                                leaf_layer = mdat_pos$leaf_layer), 
                      FUN = mean)

# invert logit back to original scale
mean_fitted_pustules <- mean(inv_log_adjusted(mdat_pos$x))  # seems plausible

mdat_pos <- mdat_pos %>% 
  dplyr::select(plot_id, position, leaf_layer, stack_id, x) %>%
  as_tibble() %>% 
  arrange(plot_id, position, leaf_layer)
# add a unique plot_layer identifier
# here, specifying each plot-layer combo as random effect seems appropriate 
# because plot and leaf_layer can be considered crossed factors, rather than nested
mdat_pos$plot_layer <- interaction(mdat_pos$plot_id, mdat_pos$leaf_layer, drop = TRUE)

# Fit AR1 along positions within plots
model_position <- lme(
  x ~ plot_id * leaf_layer,
  random = ~1 | plot_layer,
  correlation = corAR1(form = ~ as.numeric(position) | plot_layer),
  data = mdat_pos
)
plot(model_position)
qqnorm(resid(model_position))
qqline(resid(model_position))
hist(resid(model_position))
# Seems OK.

# get autocorrelation coefficient
# Valid on logit-transformed scale
phi1_pos_rust <- coef(model_position$modelStruct$corStruct, unconstrained = FALSE)

# model pieces
sigma_resid <- sqrt(model_position$sigma^2)  # residual sd
fixef_vals <- fitted(model_position, level = 0)  # fixed-effect fitted values on transformed scale
# get random effects by stack if you want to simulate them too:
ranef_df <- ranef(model_position)  # random intercepts per plot x layer

# data layout for simulation: list of stacks with their row indices
mdat_pos <- mdat_pos %>% mutate(plot_layer = fct_drop(plot_layer))
stack_groups <- split(1:nrow(mdat_pos), mdat_pos$plot_layer)

nsim <- 20
acf_sims <- matrix(NA, nrow = nsim, ncol = 25)  # lags 0..9

for(sim in 1:nsim){
  print(sim)
  sim_vals_orig <- numeric(nrow(mdat_pos))
  for(si in seq_along(stack_groups)){
    rows <- stack_groups[[si]]
    n <- length(rows)
    # Build model-predicted mean on transformed scale for these rows
    mu_trans <- fixef_vals[rows] + unname(ranef_df[as.character(mdat_pos$plot_layer[rows]), "(Intercept)"])
    # correlation matrix for this stack
    Cor <- get_cor_matrix(model_position$modelStruct$corStruct, n)
    # we multiply the correlation matrix by the residual variance to get the full covariance matrix
    # resulting sigma encodes both variance and autocorrelation for the residuals within this stack
    Sigma <- (model_position$sigma^2) * Cor
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
  acf_per_stack_mat <- do.call(rbind, acf_per_stack)  # bind vectors row-wise
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

# ============================================================================== -
# 8) get all data and make plot ----
# ============================================================================== -
model_based_acf_pos <- bind_rows(model_based_acf_pos_placl, model_based_acf_pos_rust)
model_based_acf <- bind_rows(model_based_acf_stack, model_based_acf_pos)

# add phi1
phi <- tibble(trait = c("PLACL", "PLACL", "Pustules", "Pustules"),
              level = c("stack", "pos", "stack", "pos"),
              x_pos = 0.9* c(9, 24, 9, 24),
              y_pos = c(0.9, 0.9, 0.77, 0.77),
              phi1 = c(phi1_stack_placl, phi1_pos_placl, phi1_stack_rust, phi1_pos_rust))

# add to plot
final <- both +
  geom_line(data = model_based_acf, aes(x = lag, y = mean), lty = 2) +
  # why need to add individually?
  geom_line(data = model_based_acf_pos_placl, aes(x = lag, y = mean), lty = 2) +
  geom_text(
    size = 3.5,
    data = phi,
    aes(
      x = x_pos, 
      y = y_pos,
      label = sprintf("italic(phi)[1] == %.2f", phi1)
    ),
    parse = TRUE,
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.5, 0.9))

png("all_acf_estimates.png", width = 4, height = 6, units = 'in', res = 300)
plot(final)
dev.off()

# ============================================================================== -

