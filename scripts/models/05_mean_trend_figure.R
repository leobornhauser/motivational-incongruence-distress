# Combines step 1's three per-sample mean-trend plots into one figure.
# INC is shown first, followed by GPD.
# Y-axes use the full score range for each construct.

library(ggplot2)
library(patchwork)
library(ggh4x)

source("scripts/models/01_mean_structure.R")

dir.create(
  "plots/05_mean_trend_figure",
  recursive = TRUE,
  showWarnings = FALSE
)

# Put INC before GPD in the facets
inc_first <- function(p) {
  
  facet_var <- names(p$data)[
    sapply(p$data, function(x) {
      all(c("INC", "GPD") %in% as.character(x))
    })
  ][1]
  
  p$data[[facet_var]] <- factor(
    p$data[[facet_var]],
    levels = c("INC", "GPD")
  )
  
  p
}

# Bern: INC = 1–5, GPD = 0–4
bern <- inc_first(
  plot_mean_trend(mean_trends, "Bern")
) +
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(1, 5), breaks = 1:5),
      scale_y_continuous(limits = c(0, 4), breaks = 0:4)
    )
  ) +
  labs(title = "Bern Trial: wave means (95% CI)")

# Zurich: INC = 1–5, GPD = 0–4
zurich <- inc_first(
  plot_mean_trend(mean_trends, "Zurich")
) +
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(1, 5), breaks = 1:5),
      scale_y_continuous(limits = c(0, 4), breaks = 0:4)
    )
  ) +
  labs(title = "Zurich Trial: wave means (95% CI)")

# Trier: INC = 1–5, GPD = 1–4
trier <- inc_first(
  plot_mean_trend(mean_trends, "Trier")
) +
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(1, 5), breaks = 1:5),
      scale_y_continuous(limits = c(1, 4), breaks = 1:4)
    )
  ) +
  labs(title = "Trier sample: wave means (95% CI)")

# Combine
combined <- wrap_plots(
  bern,
  zurich,
  trier,
  ncol = 1
) +
  plot_annotation(tag_levels = "A")

# Save
ggsave(
  "plots/05_mean_trend_figure/mean_trend_combined_new.png",
  plot = combined,
  width = 9,
  height = 9,
  dpi = 300,
  bg = "white"
)