library(dplyr)
library(ggplot2)
library(ggpattern)

plot_occupancy <- function(Data) {
  occ_plot_data <- Data %>%
    mutate(mined = factor(mined, levels = c("no", "yes"))) %>%
    group_by(site, spp, mined) %>%
    summarise(Z_si = as.integer(sum(count) > 0), .groups = "drop") %>%
    group_by(spp, mined) %>%
    summarise(psi_hat = mean(Z_si), .groups = "drop")

  ggplot(occ_plot_data, aes(x = spp, y = psi_hat, pattern = mined)) +
    geom_col_pattern(
      position = position_dodge(width = 0.8),
      width = 0.7,
      fill = "white",
      color = "black",
      pattern_fill = "black",
      pattern_colour = "black",
      pattern_density = 0.6,
      pattern_spacing = 0.01,
      pattern_size = 0.03
    ) +
    scale_pattern_manual(
      name = "Mining status",
      values = c("no" = "circle", "yes" = "stripe"),
      labels = c("Non-mined", "Mined")
    ) +
    labs(x = "Species", y = "Empirical occupancy probability") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

plot_abundance_heatmap <- function(Data) {
  df <- Data %>%
    group_by(site, spp, mined) %>%
    summarise(total_count = sum(count), .groups = "drop")

  site_order <- df %>%
    distinct(site, mined) %>%
    arrange(factor(mined, levels = c("no", "yes")), site) %>%
    pull(site)

  df <- df %>% mutate(site = factor(site, levels = site_order))

  site_info <- df %>%
    distinct(site, mined) %>%
    mutate(site_num = as.numeric(site))

  split_idx <- site_info %>%
    filter(mined == "yes") %>%
    summarise(min_idx = min(site_num)) %>%
    pull(min_idx)

  label_df <- site_info %>%
    group_by(mined) %>%
    summarise(x = mean(site_num), .groups = "drop") %>%
    mutate(
      label = if_else(mined == "yes", "Mined", "Non-mined"),
      y = length(unique(df$spp)) + 1
    )

  ggplot(df, aes(x = site, y = spp, fill = total_count)) +
    geom_tile() +
    geom_vline(xintercept = split_idx - 0.5, linewidth = 0.8) +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4.5
    ) +
    scale_fill_gradient(low = "white", high = "black") +
    coord_cartesian(clip = "off") +
    labs(x = "Site", y = "Species", fill = "Count") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = margin(40, 10, 10, 10)
    )
}

# Example:
# plot_occupancy(Data)
# plot_abundance_heatmap(Data)
