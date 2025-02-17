# ========================================================================= #
# PYI Impact Study                                                          #
# Plot TE SMD summary                                                       #
# Author: David Taylor                                                      #
# Date: 01/2025                                                             #
# ========================================================================= #

# This code produces a summary plot showing the results of each outcome and subgroup analysis ce intervals

# load required packages 
library(tidyverse)
library(ggplot2)
library(cowplot)

# load data
smd_plot_file_location <- "./manuscript/inputs/data/smd_plot_data.RDS"
smd_plot_data <- readRDS(smd_plot_file_location)

# declare plotting helper function to make gradient geom_rect annotations
annotate_gradient_rects <- function(n_rects, x_min, x_max, fill_color, 
                                    alpha_max = 0.2, alpha_min = 0.05,
                                    n_steps = 2, alpha_pattern = "center") {
  x_steps <- seq(from = x_min, to = x_max, length.out = n_rects + 1)
  
  # create alpha steps based on the specified pattern
  if (alpha_pattern == "center") {
    alpha_steps <- c(
      seq(from = alpha_max, to = alpha_min, length.out = n_steps),
      seq(from = alpha_min, to = alpha_max, length.out = n_steps)
    )
  } else if (alpha_pattern == "right") {  # changed from 'top'
    alpha_steps <- seq(from = alpha_max, to = alpha_min, length.out = n_steps * 2)
  } else if (alpha_pattern == "left") {   # changed from 'bottom'
    alpha_steps <- seq(from = alpha_min, to = alpha_max, length.out = n_steps * 2)
  } else {
    stop("Invalid alpha_pattern. Choose 'center', 'right', or 'left'.")
  }
  
  # repeat the alpha steps to match the number of rectangles
  alpha_steps <- rep(alpha_steps, length.out = n_rects)
  
  rect_grad <- data.frame(
    xmin = x_steps[-(n_rects + 1)],
    xmax = x_steps[-1],
    alpha = alpha_steps
  )
  
  lapply(1:nrow(rect_grad), function(i) {
    annotate("rect", 
             xmin = rect_grad$xmin[i], xmax = rect_grad$xmax[i], 
             ymin = -Inf, ymax = Inf,
             fill = fill_color, alpha = rect_grad$alpha[i])
  })
}

# subset data
smd_plot_data_filtered <- smd_plot_data |>
  filter(
    !outcome == "In new homelessness spell that requires housing assistance between 18th & 19th birthday",
    !outcome == "In new or ongoing homelessness spell that requires housing assistance between 18th & 19th birthday"
  )

# subset data for top panel
top_panel_smd_facet_plot_data <- smd_plot_data_filtered |>
  mutate(
    group = case_when(
      group == "Overall" ~ "ATT",
      group == "Sex" ~ "CATT(Sex)",
      group == "Aboriginal" ~ "CATT(Aboriginal Status)",
      group == "Housing vulnerability" ~ "CATT(Prior homelessness)",
    ),
    group = factor(
      group, 
      levels = c(
        "ATT", 
        "CATT(Sex)", 
        "CATT(Aboriginal Status)", 
        "CATT(Prior homelessness)"
      ),
      ordered = TRUE
    ),
    matching_specification = case_when(
      matching_specification == "Full" ~ "Full matching",
      matching_specification == "Nearest Neighbour" ~ "Nearest Neighbour matching"
    )
  ) |>
  filter(
    outcome %in% c(
      "In homelessness spell on 18th birthday",                                                                  
      "In homelessness spell on 19th birthday",                                                                   
      "Any homelessness spell between 18th & 19th birthday",
      "New homelessness spell between 18th & 19th birthday",
      "Days in homelessness spell between 18th and 19th birthday")
  ) |>
  mutate(
    outcome = case_when(
      outcome == "Any homelessness spell between 18th & 19th birthday" ~ "New or ongoing homelessness spell between 18th & 19th birthday",
      TRUE ~ outcome
    ),
    outcome = factor(
      outcome,
      levels = c(
        "In homelessness spell on 18th birthday",                                                                  
        "In homelessness spell on 19th birthday",                                                                   
        "New or ongoing homelessness spell between 18th & 19th birthday",
        "New homelessness spell between 18th & 19th birthday",
        "Days in homelessness spell between 18th and 19th birthday"
      ),
      ordered = TRUE
    )
  )

# subset data for bottom panel
bottom_panel_smd_facet_plot_data <- smd_plot_data_filtered |>
  mutate(
    group = case_when(
      group == "Overall" ~ "ATT",
      group == "Sex" ~ "CATT(Sex)",
      group == "Aboriginal" ~ "CATT(Aboriginal Status)",
      group == "Housing vulnerability" ~ "CATT(Prior homelessness)",
    ),
    group = factor(
      group, 
      levels = c(
        "ATT", 
        "CATT(Sex)", 
        "CATT(Aboriginal Status)", 
        "CATT(Prior homelessness)"
      ),
      ordered = TRUE
    ),
    matching_specification = case_when(
      matching_specification == "Full" ~ "Full matching",
      matching_specification == "Nearest Neighbour" ~ "Nearest Neighbour matching"
    )
  ) |>
  filter(
    outcome %in% c(
      "New unsheltered homelessness spell between 18th and 19th birthday",
      "New or ongoing unsheltered homelessness spell between 18th and 19th birthday",
      "In new homelessness spell that requires short term accommodation between 18th & 19th birthday",           
      "In new or ongoing homelessness spell that requires short term accommodation between 18th & 19th birthday",
      "Number of distinct homelessness spells between 18th and 19th birthday")
  ) |>
  mutate(
    outcome = factor(
      outcome,
      levels = c(
        "New unsheltered homelessness spell between 18th and 19th birthday",
        "New or ongoing unsheltered homelessness spell between 18th and 19th birthday",
        "In new homelessness spell that requires short term accommodation between 18th & 19th birthday",           
        "In new or ongoing homelessness spell that requires short term accommodation between 18th & 19th birthday",
        "Number of distinct homelessness spells between 18th and 19th birthday"
      ),
      ordered = TRUE
    )
  )

smd_plot_top_panel <- top_panel_smd_facet_plot_data |>  
  ggplot() +
  aes(
    x = smd,
    y = strata,
    colour = matching_specification
  ) +
  annotate_gradient_rects(
    n_rects = 50, 
    x_min = 0,
    x_max = 1,
    fill_color = "#8B0000",
    alpha_min = 0,  
    alpha_max = 0.5,  
    n_steps = 50,  
    alpha_pattern = "left"
  ) +
  annotate_gradient_rects(
    n_rects = 50, 
    x_min = 0,
    x_max = -1,
    fill_color = "#006400",
    alpha_min = 0,  
    alpha_max = 0.5,  
    n_steps = 50,  
    alpha_pattern = "left"
  ) +
  geom_vline(
    xintercept = 0,
    colour = "gray50",
    linetype = 2,
    size = 0.3
  ) +
  geom_point(
    data = top_panel_smd_facet_plot_data |> filter(matching_specification == "Full matching"),
    size = 1,
    position = position_nudge(y = 0.2)
  ) +
  geom_point(
    data = top_panel_smd_facet_plot_data |> filter(matching_specification == "Nearest Neighbour matching"),
    size = 1,
    position = position_nudge(y = -0.2)
  ) +
  geom_errorbarh(
    data = top_panel_smd_facet_plot_data |> filter(matching_specification == "Full matching"),
    position = position_nudge(y = 0.2),
    aes(
      xmin = pmax(-1, ci_low),
      xmax = pmin(1, ci_high),
      height = 0),
    linewidth = 0.25
  ) +
  geom_errorbarh(
    data = top_panel_smd_facet_plot_data |> filter(matching_specification == "Nearest Neighbour matching"),
    position = position_nudge(y = -0.2),
    aes(
      xmin = pmax(-1, ci_low),
      xmax = pmin(1, ci_high),
      height = 0),
    linewidth = 0.25
  ) +
  geom_segment(
    data = top_panel_smd_facet_plot_data |> filter(ci_low < -1 | ci_high > 1, matching_specification == "Full matching"),
    position = position_nudge(y = 0.2),
    aes(
      x = case_when(
        ci_low < -1 ~ -0.9,
        ci_high > 1 ~ 0.9
      ),
      xend = case_when(
        ci_low < -1 ~ -1,
        ci_high > 1 ~ 1
      ),
      y = strata,
      yend = strata),
    arrow = arrow(
      length = unit(0.125, "cm")),
    linewidth = 0.25
  ) +
  geom_segment(
    data = top_panel_smd_facet_plot_data |> filter(ci_low < -1 | ci_high > 1, matching_specification == "Nearest Neighbour matching"),
    position = position_nudge(y = -0.2),
    aes(
      x = case_when(
        ci_low < -1 ~ -0.9,
        ci_high > 1 ~ 0.9
      ),
      xend = case_when(
        ci_low < -1 ~ -1,
        ci_high > 1 ~ 1
      ),
      y = strata,
      yend = strata),
    arrow = arrow(
      length = unit(0.125, "cm")),
    linewidth = 0.25
  ) +
  facet_grid(
    group ~ outcome,
    scales = "free",
    space = "free",
    labeller = labeller(
      outcome = label_wrap_gen(width = 20),
      group = label_wrap_gen(width = 10)
    )
  ) +
  coord_cartesian(
    xlim = c(-1,1)
  ) +
  scale_x_continuous(
    breaks = c(-1, 0, 1)
  ) +
  scale_y_discrete(
    labels = scales::label_wrap(10)
  ) +
  scale_color_manual(
    values = c(
      "Full matching" = "#005377",
      "Nearest Neighbour matching" = "#06A77D"
    )
  ) +
  labs(
    x = "",
    y = "",
    caption = "",
    colour = "Matching specification"
  ) +
  theme(
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    strip.text.x = element_text(
      size = rel(0.7),
      colour = "#000000"
    ),
    strip.text.y = element_text(
      size = rel(0.7),
      colour = "#000000",
      face = "italic"
    ),
    axis.text.y = element_text(
      size = rel(0.7),
      colour = "#000000",
      family = "roboto"
    ),
    axis.text.x = element_text(
      size = rel(0.7),
      colour = "#000000",
      family = "roboto"
    ),
    axis.title.x = element_text(
      size = rel(0.7),
      colour = "#000000",
      family = "roboto"
    ),
    axis.title = element_text(
      size = rel(0.8),
      colour = "#000000",
      family = "roboto"
    ),
    panel.background = element_rect(
      fill = "#FFFFFF"
    ),
    legend.title = element_text(
      size = rel(0.7),
      colour = "#000000",
      family = "roboto"
    ),
    legend.text = element_text(
      size = rel(0.6),
      colour = "#000000",
      family = "roboto"
    ),
    legend.margin = margin(
      t = 0, 
      b = 0,
      unit = "pt"),
    plot.margin = margin(
      l = 0, 
      r = 0, 
      t = 0, 
      b = 0, 
      unit = "pt")
  )

smd_plot_bottom_panel <- bottom_panel_smd_facet_plot_data |>  
  ggplot() +
  aes(
    x = smd,
    y = strata,
    colour = matching_specification
  ) +
  annotate_gradient_rects(
    n_rects = 50, 
    x_min = 0,
    x_max = 1,
    fill_color = "#8B0000",
    alpha_min = 0,  
    alpha_max = 0.5,  
    n_steps = 50,  
    alpha_pattern = "left"
  ) +
  annotate_gradient_rects(
    n_rects = 50, 
    x_min = 0,
    x_max = -1,
    fill_color = "#006400",
    alpha_min = 0,  
    alpha_max = 0.5,  
    n_steps = 50,  
    alpha_pattern = "left"
  ) +
  geom_vline(
    xintercept = 0,
    colour = "gray50",
    linetype = 2,
    size = 0.3
  ) +
  geom_point(
    data = bottom_panel_smd_facet_plot_data |> filter(matching_specification == "Full matching"),
    size = 1,
    position = position_nudge(y = 0.2)
  ) +
  geom_point(
    data = bottom_panel_smd_facet_plot_data |> filter(matching_specification == "Nearest Neighbour matching"),
    size = 1,
    position = position_nudge(y = -0.2)
  ) +
  geom_errorbarh(
    data = bottom_panel_smd_facet_plot_data |> filter(matching_specification == "Full matching"),
    position = position_nudge(y = 0.2),
    aes(
      xmin = pmax(-1, ci_low),
      xmax = pmin(1, ci_high),
      height = 0),
    linewidth = 0.25
  ) +
  geom_errorbarh(
    data = bottom_panel_smd_facet_plot_data |> filter(matching_specification == "Nearest Neighbour matching"),
    position = position_nudge(y = -0.2),
    aes(
      xmin = pmax(-1, ci_low),
      xmax = pmin(1, ci_high),
      height = 0),
    linewidth = 0.25
  ) +
  geom_segment(
    data = bottom_panel_smd_facet_plot_data |> filter(ci_low < -1 | ci_high > 1, matching_specification == "Full matching"),
    position = position_nudge(y = 0.2),
    aes(
      x = case_when(
        ci_low < -1 ~ -0.9,
        ci_high > 1 ~ 0.9
      ),
      xend = case_when(
        ci_low < -1 ~ -1,
        ci_high > 1 ~ 1
      ),
      y = strata,
      yend = strata),
    arrow = arrow(
      length = unit(0.15, "cm")),
    linewidth = 0.25
  ) +
  geom_segment(
    data = bottom_panel_smd_facet_plot_data |> filter(ci_low < -1 | ci_high > 1, matching_specification == "Nearest Neighbour matching"),
    position = position_nudge(y = -0.2),
    aes(
      x = case_when(
        ci_low < -1 ~ -0.9,
        ci_high > 1 ~ 0.9
      ),
      xend = case_when(
        ci_low < -1 ~ -1,
        ci_high > 1 ~ 1
      ),
      y = strata,
      yend = strata),
    arrow = arrow(
      length = unit(0.15, "cm")),
    linewidth = 0.25
  ) +
  facet_grid(
    group ~ outcome,
    scales = "free",
    space = "free",
    labeller = labeller(
      outcome = label_wrap_gen(width = 20),
      group = label_wrap_gen(width = 10)
    )
  ) +
  coord_cartesian(
    xlim = c(-1,1)
  ) +
  scale_x_continuous(
    breaks = c(-1, 0, 1)
  ) +
  scale_y_discrete(
    labels = scales::label_wrap(10)
  ) +
  scale_color_manual(
    values = c(
      "Full matching" = "#005377",
      "Nearest Neighbour matching" = "#06A77D"
    )
  ) +
  labs(
    x = "Standardised Mean Difference",
    y = "",
    caption = "Error bars denote 95% confidence intervals",
    colour = ""
  ) +
  theme(
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    strip.text.x = element_text(
      size = rel(0.7),
      colour = "#000000"
    ),
    strip.text.y = element_text(
      size = rel(0.7),
      colour = "#000000",
      face = "italic"
    ),
    axis.text.y = element_text(
      size = rel(0.7),
      colour = "#000000",
      family = "roboto"
    ),
    axis.text.x = element_text(
      size = rel(0.7),
      colour = "#000000",
      family = "roboto"
    ),
    axis.title.x = element_text(
      size = rel(0.7),
      colour = "#000000",
      family = "roboto"
    ),
    plot.caption = element_text(
      size = rel(0.6),
      colour = "#000000",
      family = "roboto"
    ),
    panel.background = element_rect(
      fill = "#FFFFFF"
    ),
    legend.margin = margin(
      t = 0, 
      b = 0,
      unit = "pt"),
    plot.margin = margin(
      l = 0, 
      r = 0, 
      t = 0, 
      b = 0, 
      unit = "pt")
  )

# combine both plots together using cowplot
combined_smd_summary_plot <- plot_grid(
  smd_plot_top_panel, 
  smd_plot_bottom_panel, 
  ncol = 1,
  align = 'v',
  axis = 'lr',
  rel_heights = c(1, 1)
)

# add annotations
annotated_smd_summary_plot <- ggdraw() + 
  draw_plot(
    combined_smd_summary_plot
  ) +
  draw_text(
    "Favours intervention",
    x = 0.225,
    y = 0.52,
    size = 8,
    family = "roboto",
    fontface = "bold",
    color = "#006400"
  ) +
  draw_line(
    x = c(0.275, 0.15),  
    y = c(0.505, 0.505),  
    color = "#006400",
    size = 0.5,
    arrow = arrow(length = unit(0.2, "cm"))
  ) +
  draw_text(
    "Favours comparison",
    x = 0.85,
    y = 0.52,
    size = 8,
    family = "roboto",
    fontface = "bold",
    color = "#8B0000"
  ) +
  draw_line(
    x = c(0.8, 0.925),  
    y = c(0.505, 0.505),  
    color = "#8B0000",
    size = 0.5,
    arrow = arrow(length = unit(0.2, "cm"))
  )

# save plot 
ggsave(
  filename = "./analysis/output/pyi_smd_summary_plot.png",
  plot = annotated_smd_summary_plot,
  width = 6,
  height = 7.5,
  bg = "white" 
)