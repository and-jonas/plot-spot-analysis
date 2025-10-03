
library(tidyverse)
library(nlme)

setwd("O:/Data-Work/22_Plant_Production-CH/224_Digitalisation/Jonas_Anderegg_Files/C_Manuscripts/2025_FocalStack")

# ============================================================================== -
# pre-process data ----
# ============================================================================== -

d <- read_csv("overall_canopy_results.csv")

# add leaf angle to design info
ang <- c(rep("L", 6), rep("H", 6))
gen <- c("GULLIVER", "KOBRA PLUS", "HISTORY", "RAISON", "BRANDO", "POTENZIAL",
         "KWS SCIROCCO LP 509.3.04", "XENOS", "TAIFUN", "DIAVEL", "RETRO", "BORNEO")
info <- tibble(ang, gen)

# get experimental design
design <- d %>% dplyr::filter(date %in% c("2023-06-06", "2023-06-07")) %>% 
  full_join(., info, by = c("genotype_name" = "gen")) %>% 
  dplyr::select(plot_id, ang, genotype_name) %>% unique()

# select data corresponding to the largest campaign
sub <- d %>% dplyr::filter(date %in% c("2023-06-06", "2023-06-07")) %>% 
  dplyr::select(placl, `rust_density_1e-6`, plot_id, leaf_layer, position, stack_image_id)

# PLACL per genotype
pdat <- sub %>% full_join(., design)
ggplot(pdat) +
  geom_boxplot(aes(x = genotype_name, y = placl))

# # check if unique stack selected
# unique_check <- sub %>%
#   group_by(plot_id, leaf_layer, position) %>% 
#   nest() %>% 
#   mutate(n = purrr::map_dbl(data, nrow)) %>% ungroup() %>% 
#   dplyr::select(n) %>% unique()
# # does not select a unique stack for all plot x position
# 
# # try with camera id
# check_unique2 <- d %>% dplyr::filter(date %in% c("2023-06-06", "2023-06-07")) %>% 
#   mutate(cam_id = str_sub(imagefile_uid, 1, 3)) %>% 
#   dplyr::select(placl, plot_id, leaf_layer, position, stack_image_id, cam_id) %>% 
#   group_by(plot_id, leaf_layer, position, cam_id) %>% 
#   nest() %>% 
#   mutate(n = purrr::map_dbl(data, nrow))  %>% 
#   filter(n != 10)
# # does not select a unique stack for all plot x position
# 
# # remove incomplete stacks
# check_unique2 <- d %>% dplyr::filter(date %in% c("2023-06-06", "2023-06-07")) %>% 
#   mutate(cam_id = str_sub(imagefile_uid, 1, 3)) %>% 
#   dplyr::select(placl, plot_id, leaf_layer, position, stack_image_id, cam_id) %>% 
#   group_by(plot_id, leaf_layer, position, cam_id) %>% 
#   nest() %>% 
#   mutate(n = purrr::map_dbl(data, nrow))  %>% 
#   filter(n %in% c(10, 20))

# split stacks where multiple exist per plot
split_stacks <- sub %>% 
  group_by(plot_id, leaf_layer, position) %>% 
  nest() %>%
  mutate(n = purrr::map_dbl(data, nrow))  %>% 
  # remove incomplete stacks
  filter(n %in% c(10, 20)) %>% 
  mutate(data = map(data, ~ {
    n <- nrow(.x)
    if (n == 10) {
      list(.x)  # keep as single list element
    } else if (n == 20) {
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

# ============================================================================== -
# prelim: stackpos vs. placl ----
# ============================================================================== -

# visualize position-in-stack - placl relationship
# position-in-stack is associated with placl!
# positions in the beginning of the stack tend to have lower placl values than positions further down the stack
mdat <- split_stacks %>% full_join(., design) %>%
  dplyr::select(placl, `rust_density_1e-6`, plot_id, leaf_layer, position, stack_image_id, genotype_name) %>% 
  pivot_longer(placl:`rust_density_1e-6`)
pos_eff <- ggplot(mdat) +
  geom_density(aes(x = value, color = as.factor(stack_image_id)), alpha = 0.4) +
  theme(panel.grid.minor = element_blank())+
  facet_wrap(~name, scales = "free_y")
png("stackpos_traits_all.png", width = 6, height = 6, units = 'in', res = 400)
plot(pos_eff)
dev.off()

# plot for each layer
# in absolute terms, this is particularly pronounced for leaf layer 3, less for 2, and even less for 1
pdat <- sub %>% 
  filter(!(leaf_layer == 1 & placl > 0.03)) %>% 
  filter(!(leaf_layer == 2 & placl > 0.25)) %>%
  filter(!(`rust_density_1e-6` > 5)) %>%
  rename("Pustules" = 2) %>% 
  rename("PLACL" = 1) %>% 
  pivot_longer(PLACL:Pustules)

pdat <- pdat %>%
  mutate(facet_label = paste0("Layer ", leaf_layer, " - ", name))

z <- ggplot(pdat) + 
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
      nrow = 2
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = c(0.7, 0.9),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  ) +
  labs(
    x = "PLACL",   
    y = "Density"
  )

png("stackpos_traits_layers.png", width = 4, height = 6, units = 'in', res = 400)
plot(z)
dev.off()

# for each plot separately
pos_eff <- ggplot(mdat) +
  geom_density(aes(x = placl, color = as.factor(stack_image_id)), alpha = 0.4) +
  facet_wrap(~plot_id) +
  theme(panel.grid.minor = element_blank())
png("stackpos_placl_plots.png", width = 16, height = 16, units = 'in', res = 400)
plot(pos_eff)
dev.off()

# plot for each layer
pos_eff <- ggplot(mdat) +
  geom_density(aes(x = placl, color = as.factor(stack_image_id)), alpha = 0.4) +
  facet_wrap(~interaction(plot_id, leaf_layer), scales = "free") +
  theme(panel.grid.minor = element_blank())
png("stackpos_placl_plot_layers.png", width = 16, height = 32, units = 'in', res = 400)
print(pos_eff)
dev.off()

# ============================================================================== -
# acf: focal stacks ----
# ============================================================================== -

# estimate acf on un-cleaned data 
acf_list <- sub %>%
  rename("Pustules" = 2) %>% 
  rename("PLACL" = 1) %>% 
  pivot_longer(PLACL:Pustules) %>% 
  group_by(plot_id, leaf_layer, position, name) %>% 
  summarise(acf_vals = list(acf(value, plot = FALSE, lag.max = 9, na.action = na.pass)$acf),
            .groups = "drop") %>% 
  full_join(., design)

# estimate acf on cleaned data
acf_list_clean <- split_stacks %>%
  rename("Pustules" = 5) %>% 
  rename("PLACL" = 4) %>% 
  pivot_longer(PLACL:Pustules) %>% 
  group_by(plot_id, leaf_layer, position, name, stack_id) %>% 
  nest() %>% 
  mutate(n = purrr::map_dbl(data, nrow)) %>% ungroup() %>% 
  filter(n == 10) %>% unnest(cols = c(data)) %>% dplyr::select(-n) %>% 
  group_by(plot_id, leaf_layer, position, name, stack_id) %>% 
  summarise(acf_vals = list(acf(value, plot = FALSE, lag.max = 9, na.action = na.pass)$acf),
            .groups = "drop") %>% 
  full_join(., design)

# convert to long format tibble for plotting
acf_long <- acf_list_clean %>%
  mutate(acf_vals = map(acf_vals, as.numeric)) %>%
  unnest_longer(acf_vals) %>%
  group_by(plot_id, leaf_layer, position, name, stack_id) %>%
  mutate(lag = row_number() - 1) %>%
  ungroup()

# plot autocorrelation pattern for each individual stack
stacks <- ggplot(acf_long, aes(x = lag, y = acf_vals, group = interaction(plot_id, leaf_layer, position, name, stack_id))) +
  geom_line(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  labs(x = "Lag", y = "Autocorrelation", color = "Genotype", fill = "Genotype") +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  facet_wrap(~name) +
  theme(panel.grid.minor = element_blank()) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )
png("per_stack_traits.png", width = 8, height = 4, units = 'in', res = 400)
plot(stacks)
dev.off()

# # plot autocorrelation pattern for each individual stack, by leaf_layer
# stacks <- ggplot(acf_long, aes(x = lag, y = acf_vals, group = interaction(plot_id, leaf_layer, position, stack_id))) +
#   geom_line(alpha = 0.1) +
#   geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
#   labs(x = "Lag", y = "Autocorrelation", color = "Genotype", fill = "Genotype") +
#   scale_x_continuous(breaks = seq(0, 9, by = 1)) +
#   scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
#   facet_wrap(~leaf_layer) +
#   theme(panel.grid.minor = element_blank())+
#   theme(
#     strip.background = element_rect(fill = "grey90", color = NA),
#     strip.text = element_text(face = "bold"),
#     legend.position = "top",
#     panel.grid = element_blank(),
#     axis.title = element_text(face = "bold"),
#     axis.text = element_text(color = "black")
#   )
# png("per_stack_ll.png", width = 8, height = 4, units = 'in', res = 400)
# plot(stacks)
# dev.off()

# ============================================================================== -
# acf: focal stacks: mean ± SD ----
# ============================================================================== -

# get overall means and standard deviation
pdat_summary_stack <- acf_long %>%
  group_by(lag, name) %>%
  summarise(
    mean_acf = mean(acf_vals, na.rm = TRUE),
    sd_acf = sd(acf_vals, na.rm = TRUE),
    .groups = "drop"
  )

# plot mean ± SD
stack_all <- ggplot(pdat_summary_stack, aes(x = lag, y = mean_acf, group = name, color = name)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = name),
              alpha = 0.2, color = NA) +
  labs(x = "Lag", y = "Autocorrelation") +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )
png("acf_stack_all.png", width = 6, height = 6, units = 'in', res = 400)
plot(stack_all)
dev.off()
# 
# # get means and standard deviations per genotype
# pdat_summary <- acf_long %>%
#   group_by(genotype_name, lag) %>%
#   summarise(
#     mean_acf = mean(acf_vals, na.rm = TRUE),
#     sd_acf = sd(acf_vals, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # plot mean ± SD per genotype
# genos <- ggplot(pdat_summary, aes(x = lag, y = mean_acf, color = factor(genotype_name))) +
#   geom_line(size = 1) +
#   geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
#   geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = factor(genotype_name)),
#               alpha = 0.1, color = NA) +
#   labs(x = "Lag", y = "Autocorrelation", color = "Genotype", fill = "Genotype") +
#   scale_x_continuous(breaks = seq(0, 9, by = 1)) +
#   scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
#   theme(panel.grid.minor = element_blank())
# png("acf_genos.png", width = 8, height = 4, units = 'in', res = 400)
# plot(genos)
# dev.off()
# 
# # get data for genotypes
# pdat_pot <- acf_long %>%
#   group_by(genotype_name, leaf_layer, lag) %>%
#   summarise(
#     mean_acf = mean(acf_vals, na.rm = TRUE),
#     sd_acf   = sd(acf_vals, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # plot mean ± SD per leaf layer for each genotype
# geno_layer <- ggplot(pdat_pot, aes(x = lag, y = mean_acf, color = factor(leaf_layer))) +
#   geom_line(size = 1) +
#   geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
#   geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = factor(leaf_layer)),
#               alpha = 0.2, color = NA) +
#   labs(x = "Lag", y = "Autocorrelation", color = "Leaf layer", fill = "Leaf layer") +
#   scale_x_continuous(breaks = seq(0, 9, by = 1)) +
#   scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
#   facet_wrap(~genotype_name) +
#   theme(panel.grid.minor = element_blank())
# png("acf_geno_layer.png", width = 12, height = 6, units = 'in', res = 400)
# plot(geno_layer)
# dev.off()
# 
# # get means and standard deviations per leaf layer across all geontypes
# pdat_summary <- acf_long %>%
#   group_by(leaf_layer, lag) %>%
#   summarise(
#     mean_acf = mean(acf_vals, na.rm = TRUE),
#     sd_acf = sd(acf_vals, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # plot mean ± SD per leaf layer for all genotypes
# layers <- ggplot(pdat_summary, aes(x = lag, y = mean_acf, color = factor(leaf_layer))) +
#   geom_line(size = 1) +
#   geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
#   geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = factor(leaf_layer)),
#               alpha = 0.2, color = NA) +
#   labs(x = "Lag", y = "Autocorrelation", color = "Leaf layer", fill = "Leaf layer") +
#   scale_x_continuous(breaks = seq(0, 9, by = 1)) +
#   scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
#   theme(panel.grid.minor = element_blank())
# png("acf_layers.png", width = 8, height = 4, units = 'in', res = 400)
# plot(layers)
# dev.off()
# 
# # get means and standard deviations per genotype group (leaf angles)
# pdat_summary <- acf_long %>%
#   group_by(ang, lag) %>%
#   summarise(
#     mean_acf = mean(acf_vals, na.rm = TRUE),
#     sd_acf = sd(acf_vals, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # plot mean ± SD per genoype group (leaf angles)
# angles <- ggplot(pdat_summary, aes(x = lag, y = mean_acf, color = factor(ang))) +
#   geom_line(size = 1) +
#   geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
#   geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = factor(ang)),
#               alpha = 0.2, color = NA) +
#   labs(x = "Lag", y = "Autocorrelation", color = "Angle", fill = "Angle") +
#   scale_x_continuous(breaks = seq(0, 9, by = 1)) +
#   scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
#   theme(panel.grid.minor = element_blank())
# png("acf_angles.png", width = 8, height = 4, units = 'in', res = 400)
# plot(angles)
# dev.off()

# ============================================================================== -
# acf: positions ----
# ============================================================================== -

# get simple arithmetic stack means for this preliminary analysis
stack_mean <- split_stacks %>% 
  rename("Pustules" = 5) %>% 
  rename("PLACL" = 4) %>% 
  pivot_longer(PLACL:Pustules) %>% 
  group_by(plot_id, leaf_layer, position, name, stack_id) %>% 
  summarize(mean = mean(value, na.rm = T),
            .groups = "drop")

# rearrange data and estimate acf
acf_list_clean <- stack_mean %>%
  group_by(plot_id, leaf_layer, name) %>% 
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
  group_by(plot_id, leaf_layer, name) %>%
  mutate(lag = row_number() - 1) %>%
  ungroup()

# plot autocorrelation pattern for each individual stack
positions <- ggplot(acf_long, aes(x = lag, y = acf_vals, group = interaction(plot_id, leaf_layer, name))) +
  geom_line(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  labs(x = "Lag", y = "Autocorrelation", color = "Genotype", fill = "Genotype") +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  facet_wrap(~name) +
  theme(panel.grid.minor = element_blank())
png("acf_pos_plot_layer.png", width = 8, height = 4, units = 'in', res = 400)
plot(positions)
dev.off()

# ============================================================================== -
# acf: positions: mean ± SD ---- ----
# ============================================================================== -

# get overall means and standard deviation
pdat_summary_pos <- acf_long %>%
  group_by(lag, name) %>%
  summarise(
    mean_acf = mean(acf_vals, na.rm = TRUE),
    sd_acf = sd(acf_vals, na.rm = TRUE),
    .groups = "drop"
  )

# plot mean ± SD
pos_all <- ggplot(pdat_summary_pos, aes(x = lag, y = mean_acf, group = name, color = name)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = name),
              alpha = 0.2, color = NA) +
  labs(x = "Lag", y = "Autocorrelation") +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )
png("acf_pos_all.png", width = 3, height = 3, units = 'in', res = 400)
plot(pos_all)
dev.off()

# COMBINED PLOT (STACK AND POS)

pdat_summary_pos$crit <- "pos"
pdat_summary_stack$crit  <- "stack"
pdat_summary_both <- bind_rows(pdat_summary_pos, pdat_summary_stack)
both <- ggplot(pdat_summary_both, aes(x = lag, y = mean_acf, group = name, color = name)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
  geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = name),
              alpha = 0.2, color = NA) +
  labs(x = "Lag", y = "Autocorrelation") +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
  facet_wrap(~crit, ncol = 1, scales = "free_x",
             labeller = labeller(crit = c(
               "pos" = "Position in Plot",
               "stack" = "Position in Stack"
             ))) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )
png("acf_both_all.png", width = 3, height = 3, units = 'in', res = 400)
plot(both)
dev.off()

# # get means and standard deviations per leaf layer across all geontypes
# pdat_summary <- acf_long %>%
#   group_by(leaf_layer, lag) %>%
#   summarise(
#     mean_acf = mean(acf_vals, na.rm = TRUE),
#     sd_acf = sd(acf_vals, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # plot mean ± SD per leaf layer for all genotypes
# layers <- ggplot(pdat_summary, aes(x = lag, y = mean_acf, color = factor(leaf_layer))) +
#   geom_line(size = 1) +
#   geom_abline(intercept = 0, slope = 0, color = "red", lty = 2) +
#   geom_ribbon(aes(ymin = mean_acf - sd_acf, ymax = mean_acf + sd_acf, fill = factor(leaf_layer)),
#               alpha = 0.2, color = NA) +
#   labs(x = "Lag", y = "Autocorrelation", color = "Leaf layer", fill = "Leaf layer") +
#   scale_x_continuous(breaks = seq(0, 24, by = 3)) +
#   scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits = c(-0.75, 1.0)) +
#   theme(panel.grid.minor = element_blank())
# png("acf_pos_layers.png", width = 8, height = 4, units = 'in', res = 400)
# plot(layers)
# dev.off()

# ============================================================================== -
# model-based; STEP 1 ----
# Auto-correlation in focal stacks
# ============================================================================== -

# we start with a simplified scenario: data for a single position (position == 1) 
# and a single leaf layer (leaf_layer == 2) for each of the 15 plots

# subset the data
tdat <- split_stacks %>% 
  filter(stack_id == 1) %>% 
  filter(leaf_layer == 2) %>% 
  filter(position == 1) %>% 
  dplyr::select(1, 4, 6) %>% 
  # specify stack_id to identify a unique stack
  mutate(stack_id = interaction(plot_id))

# convert to factors
sub <- tdat %>%
  mutate(
    plot_id = factor(plot_id),
  )

# fit an LMM with an AR1 correlation along the positions in the stack,
# plot_id is fixed. 'stack_id' equals plot_id here. 
m_simple <- lme(
  placl ~ plot_id,
  random = ~1 | stack_id,
  correlation = corAR1(form = ~ stack_image_id | stack_id),
  data = sub,
  method = "REML",
  na.action=na.omit
)
plot(m_simple)
qqnorm(resid(m_simple))
qqline(resid(m_simple))
hist(resid(m_simple))
# VERY BAD

# logit-transformation
sub$placl_logit <- log((sub$placl + 0.001)/(1 - sub$placl + 0.001))

# re-fit the LMM
m_simple <- lme(
  placl_logit ~ plot_id, 
  random = ~1 | stack_id,
  correlation = corAR1(form = ~ stack_image_id | stack_id),
  data = sub,
  method = "REML",
  na.action=na.omit
)
plot(m_simple)
qqnorm(resid(m_simple))
qqline(resid(m_simple))
hist(resid(m_simple))
# MUCH BETTER

# ====== -

# include position (1...32) in the model. Still only one leaf layer (2).
# stack_id is now plot and position crossed.
tdat <- split_stacks %>% 
  filter(stack_id == 1) %>% 
  filter(leaf_layer == 2) %>% 
  dplyr::select(1, 3, 4, 6) %>% 
  # specify stack_id to identify a unique stack
  # Here it is the crossed plot_id and position within plot
  mutate(stack_id = interaction(plot_id, position))

# convert to factors
sub <- tdat %>%
  mutate(
    plot_id = factor(plot_id),
    position = factor(position)
  )

# logit-transformation
sub$placl_logit <- log((sub$placl + 0.001)/(1 - sub$placl + 0.001))

# Fit an LMM
# Now, we have "pseudo-replicates" from the positions (1...32) for each plot.
# Actually, we would need to specify 'position nested in plot' as a random factor to account for this. 
# However, it is apparently not possible to have different formulae for the random and the correlation term; 
# and nested terms are not allowed in the correlation term. Therefore, we include a random effect for each stack instead. 
# According to understanding, I cannot have nested random effects for plot/position and an AR(1) correlation along images in stacks at the same time. 
# So this is not perfect, but for OK for descriptive purposes. We don't do formal inference on plot-level effects ,the main goal 
# is to quantify auto-correlation based on the full data set. 
m_simple <- lme(
  placl_logit ~ plot_id, 
  random = ~1 | stack_id,
  correlation = corAR1(form = ~ stack_image_id | stack_id),
  data = sub,
  method = "REML",
  na.action=na.omit
)
plot(m_simple)
qqnorm(resid(m_simple))
qqline(resid(m_simple))
hist(resid(m_simple))

# ====== -

# include both position and leaf_layer (full data set)
tdat <- split_stacks %>% 
  filter(stack_id == 1) %>% 
  dplyr::select(1, 2, 3, 4, 6)
tdat <- tdat %>% 
  # specify stack_id to identify a unique stack
  # Here it is the crossed plot_id and position within plot
  mutate(stack_id = interaction(plot_id, position, leaf_layer))

# convert to factors
sub <- tdat %>%
  mutate(
    plot_id = factor(plot_id),
    position = factor(position),
    leaf_layer = factor(leaf_layer)
  )

# logit-transformation
sub$placl_logit <- log((sub$placl + 0.001)/(1 - sub$placl + 0.001))

# Fit an LMM
# Same issue as above. 
# Leaf layer can be considered crossed with plot (leaf layer 1 has the same/a very similar meaning across all plots)
# We include leaf layer as an additional fixed effect
m_simple0 <- lme(
  placl_logit ~ plot_id + leaf_layer,
  random = ~1 | stack_id,
  correlation = corAR1(form = ~ stack_image_id | stack_id),
  data = sub,
  method = "REML",
  na.action=na.exclude
)
plot(m_simple0)
qqnorm(resid(m_simple0))
qqline(resid(m_simple0))
hist(resid(m_simple0))

# phi1_stack <- coef(m_simple1$modelStruct$corStruct, unconstrained = FALSE)
# # AR(1) correlation for lag k: phi^k
# model_acf_stack <- tibble(
#   lag = 0:9,
#   acf_model = phi1_stack^(0:9)
# )
# model_acf_stack$crit <- "stack"


library(MASS)   # mvrnorm
library(nlme)

# helpers
inv_logit <- function(x) exp(x)/(1+exp(x))
get_cor_matrix <- function(corStruct, n){
  # corMatrix expects an initialized corStruct; we can use Initialize.corAR1
  # simpler: construct AR(1) correlation matrix from phi
  # Here handle AR(1) case:
  phi <- coef(corStruct, unconstrained = FALSE)
  mat <- phi ^ (abs(outer(1:n,1:n, "-")))
  return(mat)
}

# model pieces
sigma_resid <- sqrt(m_simple1$sigma^2)  # residual sd
fixef_vals <- fitted(m_simple1, level = 0)  # fixed-effect fitted values on transformed scale
# get random effects by stack if you want to simulate them too:
ranef_df <- ranef(m_simple1)  # random intercepts per stack

# data layout for simulation: list of stacks with their row indices
sub <- sub %>% droplevels()
stack_groups <- split(1:nrow(sub), sub$stack_id)

nsim <- 5
acf_sims <- matrix(NA, nrow = nsim, ncol = 9)  # lags 0..9

for(sim in 1:nsim){
  print(sim)
  sim_vals_orig <- numeric(nrow(sub))
  for(si in seq_along(stack_groups)){
    rows <- stack_groups[[si]]
    n <- length(rows)
    # Build model-predicted mean on transformed scale for these rows
    mu_trans <- fixef_vals[rows] + unname(ranef_df[as.character(sub$stack_id[rows]), "(Intercept)"])
    # correlation matrix for this stack
    Cor <- get_cor_matrix(m_simple1$modelStruct$corStruct, n)
    # we multiply the correlation matrix by the residual variance to get the full covariance matrix
    # resulting sigma encodes both variance and autocorrelation for the residuals within this stack
    Sigma <- (m_simple1$sigma^2) * Cor
    # simulate residuals
    eps <- mvrnorm(1, mu = rep(0,n), Sigma = Sigma)
    # add the simulated residuals to the predicted means
    sim_trans <- mu_trans + eps
    # back-transform to the original scale
    # i.e., revert the logit transformation applied to the raw data
    sim_orig <- inv_logit_adjusted(sim_trans, eps = 0.001)
    sim_vals_orig[rows] <- sim_orig
  }
  # compute mean ACF across stacks for this simulated dataset
  # compute per-stack acf and average (skip lag 0)
  acf_per_stack <- sapply(stack_groups, function(rows){
    r <- sim_vals_orig[rows]
    if(length(unique(r))>1) acf(r, lag.max=9, plot=FALSE, na.action=na.omit)$acf[-1] else rep(NA,9)
  })
  acf_per_stack_mat <- do.call(rbind, acf_per_stack)  # bind vectors row-wise
  acf_sims[sim, ] <- colMeans(acf_per_stack_mat, na.rm = TRUE)
}

# summarize
model_based_acf_mean <- colMeans(acf_sims, na.rm=TRUE)
model_based_acf_ci <- apply(acf_sims, 2, quantile, probs=c(0.025,0.975), na.rm=TRUE)








# Inverse function
inv_logit_adjusted <- function(x, eps = 0.001) {
  # Standard inverse logit
  p <- exp(x) / (1 + exp(x))
  # Undo the +eps adjustment approximately
  (p - eps) / (1 - 2*eps)
}

# Apply to residuals
resid_orig <- inv_logit_adjusted(resid(m_simple1), eps = 0.001)

# Get residuals and stack info
resid_df <- data.frame(
  stack_id = rep(names(resid(m_simple1)), 1),
  resid = as.numeric(resid(m_simple1))
)

# Make sure stack_id is a factor
resid_df$stack_id <- factor(resid_df$stack_id)

# Compute ACF for each stack separately (up to lag 9)
acf_by_stack <- resid_df %>%
  group_by(stack_id) %>%
  summarise(acf_vals = list(acf(resid, lag.max = 9, plot = FALSE)$acf))  # remove lag 0

acf_long <- acf_by_stack %>%
  unnest_longer(acf_vals) %>%
  group_by(stack_id) %>%
  mutate(lag = row_number() - 1) %>%
  ungroup()

pdat_summary_stack <- acf_long %>%
  group_by(lag) %>%
  summarise(
    mean_acf = mean(acf_vals, na.rm = TRUE),
    sd_acf = sd(acf_vals, na.rm = TRUE),
    .groups = "drop"
  )

# Compute ACF
acf(resid_orig, lag.max = 9)$acf

# Compute ACF up to lag 1
acf_res <- acf(resid_orig, lag.max = 1, plot = FALSE)

# Extract lag-1 autocorrelation
phi1_orig <- acf_res$acf[2]  # acf$acf[1] is lag 0, [2] is lag 1

phi1_orig

# ============================================================================== -
# model-based; STEP 2 ----
# Auto-correlation along positions-in-plot
# ============================================================================== -

# we obtain an estimate from the model in STEP 1 for each stack
# Get the row indices actually used by the model, subset data
used_rows <- as.numeric(rownames(m_simple1[["residuals"]]))
sub_used <- sub[used_rows, ]

# aggregate fitted values by stack
stack_means <- aggregate(fitted(m_simple1), 
                         by = list(stack_id = sub_used$stack_id,
                                   plot_id = sub_used$plot_id,
                                   position = sub_used$position,
                                   leaf_layer = sub_used$leaf_layer), 
                         FUN = mean)

# invert logit back to original scale
inv_logit <- function(x) exp(x) / (1 + exp(x))
mean_fitted_placl <- mean(inv_logit(stack_means$x))  # seems plausible

mdat <- stack_means %>% 
  dplyr::select(plot_id, position, leaf_layer, stack_id, x) %>%
  as_tibble() %>% 
  arrange(plot_id, position, leaf_layer)
mdat$plot_layer <- interaction(mdat$plot_id, mdat$leaf_layer, drop = TRUE)

# Fit AR1 along positions within plots
m_position <- lme(
  x ~ plot_id * leaf_layer,
  random = ~1 | plot_layer,
  correlation = corAR1(form = ~ as.numeric(position) | plot_layer),
  data = mdat
)
plot(m_position)
qqnorm(resid(m_position))
qqline(resid(m_position))
hist(resid(m_position))


phi1_pos <- coef(m_position$modelStruct$corStruct, unconstrained = FALSE)
# AR(1) correlation for lag k: phi^k
model_acf_pos <- tibble(
  lag = 0:24,
  acf_model = phi1_pos^(0:24)
)
model_acf_pos$crit <- "pos"
model_acf <- bind_rows(model_acf_stack, model_acf_pos)

both_with_model <- both + 
  geom_line(aes(y = acf_model, color = "Model-based"), 
            data = model_acf, size = 1) +
  geom_line(aes(y = mean_acf, color = "Sample-based"), 
            data = pdat_summary_both, size = 1) +
  scale_color_manual(values = c("Sample-based" = "darkblue", "Model-based" = "darkgreen")) +
  theme(
    legend.position = c(0.2, 0.9),  # x = 0.8, y = 0.9 in panel coordinates
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.title = element_blank()
  )

png("acf_both_all_with_model.png", width = 4, height = 6, units = 'in', res = 400)
plot(both_with_model)
dev.off()

# ============================================================================== -

# in one go: stack

# 1) PLACL

# include both position and leaf_layer (full data set)
tdat_placl <- split_stacks %>% 
  filter(stack_id == 1) %>% 
  dplyr::select(1, 2, 3, 4, 6)
tdat_placl <- tdat_placl %>% 
  # specify stack_id to identify a unique stack
  # Here it is the crossed plot_id and position within plot
  mutate(stack_id = interaction(plot_id, position, leaf_layer))

# convert to factors
sub_placl <- tdat_placl %>%
  mutate(
    plot_id = factor(plot_id),
    position = factor(position),
    leaf_layer = factor(leaf_layer)
  )

# logit-transformation
sub_placl$placl_logit <- log((sub$placl + 0.001)/(1 - sub$placl + 0.001))

m_simple1 <- lme(
  placl_logit ~ plot_id * leaf_layer,
  random = ~1 | stack_id,
  correlation = corAR1(form = ~ stack_image_id | stack_id),
  data = sub_placl,
  method = "REML",
  na.action=na.omit
)
plot(m_simple1)
qqnorm(resid(m_simple1))
qqline(resid(m_simple1))
hist(resid(m_simple1))
phi1_stack <- coef(m_simple1$modelStruct$corStruct, unconstrained = FALSE)

# 2) RUST

# include both position and leaf_layer (full data set)
tdat_rust <- split_stacks %>% 
  filter(stack_id == 1) %>% 
  dplyr::select(1, 2, 3, 5, 6)
tdat_rust <- tdat_rust%>% 
  # specify stack_id to identify a unique stack
  # Here it is the crossed plot_id and position within plot
  mutate(stack_id = interaction(plot_id, position, leaf_layer))

# convert to factors
sub_rust <- tdat_rust %>%
  mutate(
    plot_id = factor(plot_id),
    position = factor(position),
    leaf_layer = factor(leaf_layer)
  )

# logit-transformation
names(sub_rust)[4] <- "rust"
sub_rust <- sub_rust %>% 
  filter(complete.cases(.))
ggplot(sub_rust) +
  geom_histogram(aes(x = rust), bins = 100) +
  scale_x_continuous(limits = c(0, 50))
sub_rust$rust_log <- log1p(sub_rust$rust)
ggplot(sub_rust) +
  geom_histogram(aes(x = rust_log), bins = 50)

m_simple1 <- lme(
  rust_log ~ plot_id * leaf_layer,
  random = ~1 | stack_id,
  correlation = corAR1(form = ~ stack_image_id | stack_id),
  data = sub_rust,
  method = "REML",
  na.action=na.exclude
)
plot(m_simple1)
qqnorm(resid(m_simple1))
qqline(resid(m_simple1))
hist(resid(m_simple1))

m_simple1_varp <- lme(
  rust_log ~ plot_id * leaf_layer,
  random = ~1 | stack_id,
  correlation = corAR1(form = ~ stack_image_id | stack_id),
  data = sub_rust,
  method = "REML",
  weights = varPower(form = ~ fitted(.)),   # models variance increasing with fitted values
  na.action=na.exclude
)

AIC(m_simple1, m_simple1_varp)

plot(m_simple1_varp)
qqnorm(resid(m_simple1_varp))
qqline(resid(m_simple1_varp))
hist(resid(m_simple1_varp))
phi1_stack <- coef(m_simple1_varp$modelStruct$corStruct, unconstrained = FALSE)

fitted_vals <- fitted(m_simple1_varp)
resid_vals  <- resid(m_simple1_varp, type = "pearson")

plot_df <- data.frame(fitted = fitted_vals, resid = resid_vals)
ggplot(plot_df, aes(x = fitted, y = resid^2)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  labs(y = "Squared residuals", x = "Fitted values") +
  theme_classic()
