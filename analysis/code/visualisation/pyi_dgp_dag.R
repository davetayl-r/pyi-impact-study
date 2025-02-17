library(tidyverse)
library(ggdag)
library(dagitty)
library(cowplot)

# specify dag
dag_specification <- dagify(
  # specify estimand
  homeless ~ pyi,
  # factors which affect selection into the intervention
  pyi ~ time_care, 
  pyi ~ instab, 
  pyi ~ resi, 
  pyi ~ perman,
  # factors which affect the outcome
  homeless ~ self_placed,
  homeless ~ instab, 
  homeless ~ resi, 
  homeless ~ prior_shs,
  homeless ~ kinship,
  homeless ~ time_care,
  homeless ~ perman,
  homeless ~ unobserved_before_18,
  # factors which are affected by the intervention
  unobserved_post_18_mediator ~ pyi,
  unobserved_post_18_compete ~ pyi,
  # relationships between unobserved mediator/outcomes
  unobserved_post_18_compete ~ unobserved_post_18_mediator,
  # factors which affect intervention / outcome mediator
  homeless ~ unobserved_post_18_mediator,
  unobserved_post_18_mediator ~ self_placed,
  unobserved_post_18_mediator ~ instab, 
  unobserved_post_18_mediator ~ resi, 
  unobserved_post_18_mediator ~ prior_shs,
  unobserved_post_18_mediator ~ kinship,
  unobserved_post_18_mediator ~ time_care,
  unobserved_post_18_mediator ~ perman,
  unobserved_post_18_mediator ~ unobserved_before_18,
  # factors that affect competing outcomes
  unobserved_post_18_compete ~ homeless,
  unobserved_post_18_compete ~ self_placed,
  unobserved_post_18_compete ~ instab, 
  unobserved_post_18_compete ~ resi, 
  unobserved_post_18_compete ~ prior_shs,
  unobserved_post_18_compete ~ kinship,
  unobserved_post_18_compete ~ time_care,
  unobserved_post_18_compete ~ perman,
  unobserved_post_18_compete ~ unobserved_before_18,
  # specify potential relationships between confounds
  homeless ~ unobserved_oohc,
  prm ~ unobserved_oohc,
  time_care ~ unobserved_oohc,
  kinship ~ unobserved_oohc,
  instab ~ unobserved_oohc,
  resi ~ unobserved_oohc,
  perman ~ unobserved_oohc,
  self_placed ~ unobserved_oohc,
  prior_shs ~ unobserved_oohc,
  unobserved_before_18 ~ unobserved_oohc,
  kinship ~ time_care,
  instab ~ time_care,
  resi ~ time_care,
  perman ~ time_care,
  self_placed ~ time_care,
  prior_shs ~ time_care,
  unobserved_before_18 ~ time_care,
  kinship ~ instab,
  resi ~ instab,
  perman ~ instab,
  self_placed ~ instab,
  prior_shs ~ instab,
  unobserved_before_18 ~ instab,
  kinship ~ prior_shs,
  resi ~ prior_shs,
  perman ~ prior_shs,
  self_placed ~ prior_shs,
  unobserved_before_18 ~ prior_shs,
  self_placed ~ kinship,
  unobserved_before_18 ~ kinship,
  perman ~ resi,
  self_placed ~ resi,
  unobserved_before_18 ~ resi,
  self_placed ~ perman,
  unobserved_before_18 ~ self_placed,
  unobserved_before_18 ~ perman,
  time_care ~ prm,
  kinship ~ prm,
  resi ~ prm,
  perman ~ prm,
  self_placed ~ prm,
  unobserved_before_18 ~ prm,
  # factors that affect competing outcomes/mediators
  unobserved_post_18_mediator ~ unobserved_before_18,
  # assign labels  
  labels = c(
    "homeless" = "Y[1]", 
    "pyi" = "D",
    "unobserved_oohc" = "Z[1]",
    "time_care" = "X[1]", 
    "prm" = "Z[2]",
    "instab" = "X[2]",
    "resi" = "X[3]",
    "perman" = "X[4]",
    "kinship" = "Z[3]",
    "self_placed" = "Z[4]",
    "prior_shs" = "Z[5]",
    "unobserved_before_18" = "Z[6]",
    "unobserved_post_18_mediator" = "M",
    "unobserved_post_18_compete" = "Y[2]"
  ),
  # set coords
  coords = data.frame(
    name = c(
      "homeless",
      "pyi",
      "unobserved_oohc",
      "time_care",
      "prm",
      "instab",
      "prior_shs",
      "kinship",
      "resi",
      "perman",
      "self_placed",
      "unobserved_before_18",
      "unobserved_post_18_mediator",
      "unobserved_post_18_compete"
    ),
    x = c(
      15, # homeless
      11, # pyi
      1, # unobserved_oohc
      2, # time_care
      2, # prm
      3, # instab
      4, # prior_shs
      5, # kinship
      5, # resi
      6, # perman
      7, # self_placed
      8, # unobserved_before_18
      13, # unobserved_post_18_mediator
      15 # unobserved_post_18_compete
    ),
    y = c(
      0, # homeless
      0, # pyi
      9, # unobserved_oohc
      8, # time_care
      1, # prm
      7, # instab
      6, # prior_shs
      5, # kinship
      4, # resi
      3, # perman
      2, # self_placed
      1, # unobserved_before_18
      -0, # unobserved_post_18_mediator
      -2 # unobserved_post_18_compete
    )
  ),
  exposure = "pyi",
  outcome = "homeless"
)

# check if adjustment specific set is valid
dagitty::isAdjustmentSet(
  dag_specification,
  Z = c("instab", "perman", "resi", "time_care", "prior_shs", "kinship", "prm"),
  exposure = "pyi",
  outcome = "homeless" 
)

# prep plot data
pyi_dag_plot_data <- dag_specification |>
  tidy_dagitty() |>
  # name edges
  mutate(
    edge_name = paste(name, to, sep = "_")
  ) |>
  # specify node type
  mutate(
    node_type = case_when(
      str_detect(name, "unobserved") ~ "unobserved",
      TRUE ~ "observed"),
    var_type = case_when(
      name == "homeless" ~ "outcome",
      name == "pyi" ~ "intervention",
      name == "instab" ~ "adjusted",
      name == "time_care" ~ "adjusted",
      name == "prior_shs" ~ "adjusted",
      name == "kinship" ~ "adjusted",
      name == "resi" ~ "adjusted",
      name == "perman" ~ "adjusted",
      name == "prm" ~ "adjusted",
      TRUE ~ "unadjusted"
    ),
    edge_colour = case_when(
      var_type == "adjusted" ~ "#0EAD69",
      TRUE ~ "#3C3C3C"
    ),
    edge_alpha = case_when(
      var_type == "adjusted" ~ 0.5,
      TRUE ~ 0.2
  ))

# plot dag
pyi_dag <- pyi_dag_plot_data |>
  ggplot(
    aes(
      x = x, 
      y = y, 
      xend = xend, 
      yend = yend)) +
  geom_dag_edges_arc(
    aes(
      edge_alpha = edge_alpha,
      edge_colour = edge_colour
      ),
    curvature = c(
      homeless_unobserved_post_18_compete = 0,
      instab_homeless = 0.3,
      instab_kinship = -0.6,
      instab_perman = -0.6,
      instab_prior_shs = -0.6,
      instab_pyi = 0.2,
      instab_resi = -0.6,
      instab_self_placed = -0.6,
      instab_unobserved_before_18 = -0.6,
      instab_unobserved_post_18_compete = -0.5,
      instab_unobserved_post_18_mediator = 0.3,
      kinship_homeless = 0.3,
      kinship_self_placed = -0.6,
      kinship_unobserved_before_18 = -0.6,
      kinship_unobserved_post_18_compete = -0.5,
      kinship_unobserved_post_18_mediator = 0.3,
      perman_homeless = 0.3,
      perman_pyi = 0.2,
      perman_self_placed = -0.6,
      perman_unobserved_before_18 = -0.6,
      perman_unobserved_post_18_compete = -0.35,
      perman_unobserved_post_18_mediator = 0.3,
      prior_shs_homeless = 0.3,
      prior_shs_kinship = -0.6,
      prior_shs_perman = -0.6,
      prior_shs_resi = -0.6,
      prior_shs_self_placed = -0.6,
      prior_shs_unobserved_before_18 = -0.6,
      prior_shs_unobserved_post_18_compete = -0.5,
      prior_shs_unobserved_post_18_mediator = 0.3,
      prm_kinship = 0.1,
      prm_perman = -0.1,
      prm_resi = -0.1,
      prm_self_placed = -0.1,
      prm_time_care = 0.1,
      prm_unobserved_before_18 = -0.3,
      pyi_homeless = -0.5,
      pyi_unobserved_post_18_compete = -0.2,
      pyi_unobserved_post_18_mediator = 0,
      resi_homeless = 0.3,
      resi_perman = -0.6,
      resi_pyi = 0.2,
      resi_self_placed = -0.6,
      resi_unobserved_before_18 = -0.6,
      resi_unobserved_post_18_compete = -0.4,
      resi_unobserved_post_18_mediator = 0.3,
      self_placed_homeless = 0.3,
      self_placed_unobserved_before_18 = -0.6,
      self_placed_unobserved_post_18_compete = -0.3,
      self_placed_unobserved_post_18_mediator = 0.3,
      time_care_homeless = 0.3,
      time_care_instab = -0.6,
      time_care_kinship = -0.6,
      time_care_perman = -0.6,
      time_care_prior_shs = -0.6,
      time_care_pyi = 0.3,
      time_care_resi = -0.6,
      time_care_self_placed = -0.6,
      time_care_unobserved_before_18 = -0.6,
      time_care_unobserved_post_18_compete = -0.5,
      time_care_unobserved_post_18_mediator = 0.3,
      unobserved_before_18_homeless = 0.1,
      unobserved_before_18_unobserved_post_18_compete = -0.2,
      unobserved_before_18_unobserved_post_18_mediator = 0.3,
      unobserved_oohc_homeless = 0.3,
      unobserved_oohc_instab = -0.6,
      unobserved_oohc_kinship = -0.6,
      unobserved_oohc_perman = -0.6,
      unobserved_oohc_prior_shs = -0.6,
      unobserved_oohc_prm = -0.2,
      unobserved_oohc_resi = -0.6,
      unobserved_oohc_self_placed = -0.6,
      unobserved_oohc_time_care = -0.6,
      unobserved_oohc_unobserved_before_18 = -0.6,
      unobserved_post_18_compete_NA = 0,
      unobserved_post_18_mediator_homeless = 0
    )
  ) +
  geom_dag_node(
    aes(
      shape = node_type,
      colour = var_type),
    size = 12) +
  geom_dag_text(
    aes(
      label = label
    ),
    parse = TRUE
  ) +
  scale_colour_manual(
    values = c(
      adjusted = "#0EAD69", 
      outcome = "#005377",
      intervention = "#EE4266",
      unadjusted = "#3C3C3C"),
    breaks = c(
      "intervention", 
      "outcome", 
      "adjusted", 
      "unadjusted"),
    labels = c(
      "outcome" = "Outcome",
      "intervention" = "Intervention",
      "adjusted" = "Adjusted", 
      "unadjusted" = "Unadjusted"),
    name = "Node:"
  ) +
  scale_shape_manual(
    values = c(
      unobserved = 18, 
      observed = 16
      ),
    labels = c(
      "observed" = "Observed", 
      "unobserved" = "Unobserved"),
    name = "Variable:"
  ) +
  theme_void() +
  theme(
    legend.box = "vertical",  
    legend.position = "bottom",
    legend.title = element_text(
      face = "italic"
      )
  )

# create legend text plot
legend_plot <- 
  ggplot() + 
  scale_x_continuous(
    limits = c(0, 1)
    ) +
  scale_y_continuous(
    limits = c(-1.2, 0.2)
    ) +
  annotate("text",
    x = 0, 
    y = 0,
    hjust = 0,
    label = "Key:",
    fontface = "italic",
    size = rel(4)
    ) +
  annotate("text",
    x = 0, 
    y = -0.25,
    hjust = 0,
    size = rel(3),
    label = expression(paste(
      Y[1], ": Use of SHS after age 18, ",
      Y[2], ": Other (unobserved) outcomes, ",
      D, ": Intervention (PYI), ",
      M, ": Mediators",
    ))
  ) +
  annotate("text",
    x = 0, 
    y = -0.5,
    hjust = 0,
    size = rel(3),
    label = expression(paste(
      X[1], ": 12 months or more in OOHC, ",
      X[2], ": History of placement instability, ", 
      X[3], ": In residential care placement during eligibility period, "
    ))
  ) +
  annotate("text",
    x = 0, 
    y = -0.75,
    hjust = 0,
    size = rel(3),
    label = expression(paste(
      X[4], ": In permanent care placement during eligibility period, ",
      Z[1], ": Factors that occured before or during OOHC, ",
      Z[2], ": Parental responsibility of the Minister during eligibility period, "
    )) 
  ) +
  annotate("text",
    x = 0, 
    y = -1,
    hjust = 0,
    size = rel(3),
    label = expression(paste(
      Z[3], ": In kinship care placement during eligibility period, ",
      Z[4], ": Self-placed from placement after age 16, ",
      Z[5], ": Use of SHS between age 16 and 18, ",
      Z[6], ": Factors that occured before 18"
    ))
  ) +
  theme_void()

# bring together dag and legend
export_dag <- plot_grid(
  pyi_dag,
  legend_plot,
  ncol = 1,
  rel_heights = c(0.85, 0.15)
)

# save plot 
ggsave(
  filename = "./analysis/output/pyi_dgp_dag.png",
  plot = export_dag,
  width = 10,
  height = 7,
  bg = "white" 
)