
library(tidyverse)
library(viridis)

compute_mean_sem_by_group <- function(csv_file,
                                      group_cols,
                                      value_col,
                                      sem_plot = "sem_plot.png",
                                      std_plot = "std_plot.png") {
  
  # Read data
  df <- read_csv(csv_file)
  
  # Group and summarize
  summary_df <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      count = n(),
      max_stack = max(stack_image_id, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      sem = sd / sqrt(count),
      adjusted_count = count / max_stack
    )
  
  #------------------------------------------------------
  # SEM plot
  #------------------------------------------------------
  
  p_sem <- ggplot(summary_df,
                  aes(mean, sem, colour = adjusted_count)) +
    geom_point(size = 2.5) +
    scale_colour_viridis_c(trans = "log10") +
    labs(
      x = "Mean PLACL",
      y = "Standard Error of the Mean (PLACL)",
      colour = "Number of\nScanning Positions"
    ) +
    theme_bw(base_size = 14)
  
  ggsave(
    sem_plot,
    p_sem,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  #------------------------------------------------------
  # SD plot
  #------------------------------------------------------
  
  p_sd <- ggplot(summary_df,
                 aes(mean, sd, colour = adjusted_count)) +
    geom_point(size = 2.5) +
    scale_colour_viridis_c(trans = "log10") +
    labs(
      x = "Mean PLACL",
      y = "Standard Deviation of PLACL",
      colour = "Number of\nScanning Positions"
    ) +
    theme_bw(base_size = 14)
  
  ggsave(
    std_plot,
    p_sd,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  summary_df
}


get_top_adjusted_counts <- function(csv_file,
                                    group_cols,
                                    value_col) {
  
  df <- read_csv(csv_file)
  
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      count = n(),
      max_stack = max(stack_image_id, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      adjusted_count = count / max_stack
    ) %>%
    arrange(desc(adjusted_count)) %>%
    slice_head(n = 10)
}

#------------------------------------------------------
# Example
#------------------------------------------------------

results <- compute_mean_sem_by_group(
  csv_file = "/projects/paper2/main_run/overall_canopy_results.csv",
  group_cols = c("plot_id", "leaf_layer", "date"),
  value_col = "placl"
)

top10 <- get_top_adjusted_counts(
  csv_file = "/projects/paper2/main_run/overall_canopy_results.csv",
  group_cols = c("plot_id", "leaf_layer", "date"),
  value_col = "placl"
)

print(results)
print(top10)