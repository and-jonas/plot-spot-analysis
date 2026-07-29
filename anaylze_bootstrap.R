
library(tidyverse)


d <- readRDS("O:/Data-Work/22_Plant_Production-CH/224_Digitalisation/Jonas_Anderegg_Files/C_Manuscripts/2025_FocalStack_Downstream/data/bootstrap_results_placl.rds")

n_eff_foc_placl <- lapply(d, "[[", 1) %>% unlist()
n_eff_foc_placl_mean <- mean(n_eff_foc_placl)
quantile(
  n_eff_foc_placl,
  probs = c(0.025, 0.975)
)

n_eff_pos_placl <- lapply(d, "[[", 2) %>% unlist()
n_eff_pos_placl_mean <- mean(n_eff_pos_placl)
quantile(
  n_eff_pos_placl,
  probs = c(0.025, 0.975)
)

