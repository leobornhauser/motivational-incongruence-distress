# Step 3 of 4: the project's main, detrended RI-CLPM (Figure 1) -- fits
# unconstrained and constrained (equal-over-time) autoregressive/
# cross-lagged paths. 02 riclpm mean structure/ sources this script


library(readr)
library(dplyr)
library(psych)
library(tibble)
library(stringr)
library(lavaan)
library(semTools)
library(semPlot)
library(tidySEM)
library(ggplot2)

dir.create("plots/03 riclpm", recursive = TRUE, showWarnings = FALSE)

########## Data ##################################################################

load("data/data.RData")

head(df_bern)
head(df_zurich)
head(df_trier)

########## Define the RI-CLPM: 4 waves (Bern and Zurich) ########################

# Unconstrained: autoregressive and cross-lagged paths are free to differ
# across waves.
model_riclpm_4w <- '

  # random-intercept factors
  RI_inc =~ 1*ink1 + 1*ink2 + 1*ink3 + 1*ink4
  RI_gpd =~ 1*gpd1 + 1*gpd2 + 1*gpd3 + 1*gpd4

  # wave-specific residual factors
  inc_w_t1 =~ 1*ink1
  inc_w_t2 =~ 1*ink2
  inc_w_t3 =~ 1*ink3
  inc_w_t4 =~ 1*ink4

  gpd_w_t1 =~ 1*gpd1
  gpd_w_t2 =~ 1*gpd2
  gpd_w_t3 =~ 1*gpd3
  gpd_w_t4 =~ 1*gpd4

  # the observed scores have no residual variance of their own: everything
  # is carried by the random intercept and the wave-specific factor above
  ink1 ~~ 0*ink1
  ink2 ~~ 0*ink2
  ink3 ~~ 0*ink3
  ink4 ~~ 0*ink4

  gpd1 ~~ 0*gpd1
  gpd2 ~~ 0*gpd2
  gpd3 ~~ 0*gpd3
  gpd4 ~~ 0*gpd4

  # structural relations between residual factors
  # autoregressive (stability) effects
  inc_w_t2 ~ inc_w_t1
  inc_w_t3 ~ inc_w_t2
  inc_w_t4 ~ inc_w_t3

  gpd_w_t2 ~ gpd_w_t1
  gpd_w_t3 ~ gpd_w_t2
  gpd_w_t4 ~ gpd_w_t3

  # cross-lagged effects
  gpd_w_t2 ~ inc_w_t1
  gpd_w_t3 ~ inc_w_t2
  gpd_w_t4 ~ inc_w_t3

  inc_w_t2 ~ gpd_w_t1
  inc_w_t3 ~ gpd_w_t2
  inc_w_t4 ~ gpd_w_t3

  # within-person covariances at each wave
  inc_w_t1 ~~ gpd_w_t1
  inc_w_t2 ~~ gpd_w_t2
  inc_w_t3 ~~ gpd_w_t3
  inc_w_t4 ~~ gpd_w_t4

  # random-intercept factors: they may covary with each other, but by
  # construction are unrelated to the wave-one residual factors
  RI_inc ~~ RI_gpd
  RI_inc ~~ 0*inc_w_t1
  RI_inc ~~ 0*gpd_w_t1
  RI_gpd ~~ 0*inc_w_t1
  RI_gpd ~~ 0*gpd_w_t1
'

# Constrained: autoregressive (labels a, b) and cross-lagged (labels c, d)
# paths held equal across waves. This is the model in Figure 1.
model_riclpm_4w_constrained <- '

  # random-intercept factors
  RI_inc =~ 1*ink1 + 1*ink2 + 1*ink3 + 1*ink4
  RI_gpd =~ 1*gpd1 + 1*gpd2 + 1*gpd3 + 1*gpd4

  # wave-specific residual factors
  inc_w_t1 =~ 1*ink1
  inc_w_t2 =~ 1*ink2
  inc_w_t3 =~ 1*ink3
  inc_w_t4 =~ 1*ink4

  gpd_w_t1 =~ 1*gpd1
  gpd_w_t2 =~ 1*gpd2
  gpd_w_t3 =~ 1*gpd3
  gpd_w_t4 =~ 1*gpd4

  ink1 ~~ 0*ink1
  ink2 ~~ 0*ink2
  ink3 ~~ 0*ink3
  ink4 ~~ 0*ink4

  gpd1 ~~ 0*gpd1
  gpd2 ~~ 0*gpd2
  gpd3 ~~ 0*gpd3
  gpd4 ~~ 0*gpd4

  # structural relations between residual factors
  # autoregressive effects, held equal across waves
  inc_w_t2 ~ a*inc_w_t1
  inc_w_t3 ~ a*inc_w_t2
  inc_w_t4 ~ a*inc_w_t3

  gpd_w_t2 ~ b*gpd_w_t1
  gpd_w_t3 ~ b*gpd_w_t2
  gpd_w_t4 ~ b*gpd_w_t3

  # cross-lagged effects, held equal across waves
  gpd_w_t2 ~ c*inc_w_t1
  gpd_w_t3 ~ c*inc_w_t2
  gpd_w_t4 ~ c*inc_w_t3

  inc_w_t2 ~ d*gpd_w_t1
  inc_w_t3 ~ d*gpd_w_t2
  inc_w_t4 ~ d*gpd_w_t3

  inc_w_t1 ~~ gpd_w_t1
  inc_w_t2 ~~ gpd_w_t2
  inc_w_t3 ~~ gpd_w_t3
  inc_w_t4 ~~ gpd_w_t4

  # random-intercept factors: they may covary with each other, but by
  # construction are unrelated to the wave-one residual factors
  RI_inc ~~ RI_gpd
  RI_inc ~~ 0*inc_w_t1
  RI_inc ~~ 0*gpd_w_t1
  RI_gpd ~~ 0*inc_w_t1
  RI_gpd ~~ 0*gpd_w_t1
'

########## Define the RI-CLPM: 6 waves (Trier) ###################################

model_riclpm_6w <- '

  # random-intercept factors
  RI_inc =~ 1*ink1 + 1*ink2 + 1*ink3 + 1*ink4 + 1*ink5 + 1*ink6
  RI_gpd =~ 1*gpd1 + 1*gpd2 + 1*gpd3 + 1*gpd4 + 1*gpd5 + 1*gpd6

  # wave-specific residual factors
  inc_w_t1 =~ 1*ink1
  inc_w_t2 =~ 1*ink2
  inc_w_t3 =~ 1*ink3
  inc_w_t4 =~ 1*ink4
  inc_w_t5 =~ 1*ink5
  inc_w_t6 =~ 1*ink6

  gpd_w_t1 =~ 1*gpd1
  gpd_w_t2 =~ 1*gpd2
  gpd_w_t3 =~ 1*gpd3
  gpd_w_t4 =~ 1*gpd4
  gpd_w_t5 =~ 1*gpd5
  gpd_w_t6 =~ 1*gpd6

  ink1 ~~ 0*ink1
  ink2 ~~ 0*ink2
  ink3 ~~ 0*ink3
  ink4 ~~ 0*ink4
  ink5 ~~ 0*ink5
  ink6 ~~ 0*ink6

  gpd1 ~~ 0*gpd1
  gpd2 ~~ 0*gpd2
  gpd3 ~~ 0*gpd3
  gpd4 ~~ 0*gpd4
  gpd5 ~~ 0*gpd5
  gpd6 ~~ 0*gpd6

  # structural relations between residual factors
  inc_w_t2 ~ inc_w_t1
  inc_w_t3 ~ inc_w_t2
  inc_w_t4 ~ inc_w_t3
  inc_w_t5 ~ inc_w_t4
  inc_w_t6 ~ inc_w_t5

  gpd_w_t2 ~ gpd_w_t1
  gpd_w_t3 ~ gpd_w_t2
  gpd_w_t4 ~ gpd_w_t3
  gpd_w_t5 ~ gpd_w_t4
  gpd_w_t6 ~ gpd_w_t5

  gpd_w_t2 ~ inc_w_t1
  gpd_w_t3 ~ inc_w_t2
  gpd_w_t4 ~ inc_w_t3
  gpd_w_t5 ~ inc_w_t4
  gpd_w_t6 ~ inc_w_t5

  inc_w_t2 ~ gpd_w_t1
  inc_w_t3 ~ gpd_w_t2
  inc_w_t4 ~ gpd_w_t3
  inc_w_t5 ~ gpd_w_t4
  inc_w_t6 ~ gpd_w_t5

  inc_w_t1 ~~ gpd_w_t1
  inc_w_t2 ~~ gpd_w_t2
  inc_w_t3 ~~ gpd_w_t3
  inc_w_t4 ~~ gpd_w_t4
  inc_w_t5 ~~ gpd_w_t5
  inc_w_t6 ~~ gpd_w_t6

  # random-intercept factors: they may covary with each other, but by
  # construction are unrelated to the wave-one residual factors
  RI_inc ~~ RI_gpd
  RI_inc ~~ 0*inc_w_t1
  RI_inc ~~ 0*gpd_w_t1
  RI_gpd ~~ 0*inc_w_t1
  RI_gpd ~~ 0*gpd_w_t1
'

model_riclpm_6w_constrained <- '

  # random-intercept factors
  RI_inc =~ 1*ink1 + 1*ink2 + 1*ink3 + 1*ink4 + 1*ink5 + 1*ink6
  RI_gpd =~ 1*gpd1 + 1*gpd2 + 1*gpd3 + 1*gpd4 + 1*gpd5 + 1*gpd6

  # wave-specific residual factors
  inc_w_t1 =~ 1*ink1
  inc_w_t2 =~ 1*ink2
  inc_w_t3 =~ 1*ink3
  inc_w_t4 =~ 1*ink4
  inc_w_t5 =~ 1*ink5
  inc_w_t6 =~ 1*ink6

  gpd_w_t1 =~ 1*gpd1
  gpd_w_t2 =~ 1*gpd2
  gpd_w_t3 =~ 1*gpd3
  gpd_w_t4 =~ 1*gpd4
  gpd_w_t5 =~ 1*gpd5
  gpd_w_t6 =~ 1*gpd6

  ink1 ~~ 0*ink1
  ink2 ~~ 0*ink2
  ink3 ~~ 0*ink3
  ink4 ~~ 0*ink4
  ink5 ~~ 0*ink5
  ink6 ~~ 0*ink6

  gpd1 ~~ 0*gpd1
  gpd2 ~~ 0*gpd2
  gpd3 ~~ 0*gpd3
  gpd4 ~~ 0*gpd4
  gpd5 ~~ 0*gpd5
  gpd6 ~~ 0*gpd6

  # structural relations between residual factors
  inc_w_t2 ~ a*inc_w_t1
  inc_w_t3 ~ a*inc_w_t2
  inc_w_t4 ~ a*inc_w_t3
  inc_w_t5 ~ a*inc_w_t4
  inc_w_t6 ~ a*inc_w_t5

  gpd_w_t2 ~ b*gpd_w_t1
  gpd_w_t3 ~ b*gpd_w_t2
  gpd_w_t4 ~ b*gpd_w_t3
  gpd_w_t5 ~ b*gpd_w_t4
  gpd_w_t6 ~ b*gpd_w_t5

  gpd_w_t2 ~ c*inc_w_t1
  gpd_w_t3 ~ c*inc_w_t2
  gpd_w_t4 ~ c*inc_w_t3
  gpd_w_t5 ~ c*inc_w_t4
  gpd_w_t6 ~ c*inc_w_t5

  inc_w_t2 ~ d*gpd_w_t1
  inc_w_t3 ~ d*gpd_w_t2
  inc_w_t4 ~ d*gpd_w_t3
  inc_w_t5 ~ d*gpd_w_t4
  inc_w_t6 ~ d*gpd_w_t5

  inc_w_t1 ~~ gpd_w_t1
  inc_w_t2 ~~ gpd_w_t2
  inc_w_t3 ~~ gpd_w_t3
  inc_w_t4 ~~ gpd_w_t4
  inc_w_t5 ~~ gpd_w_t5
  inc_w_t6 ~~ gpd_w_t6

  # random-intercept factors: they may covary with each other, but by
  # construction are unrelated to the wave-one residual factors
  RI_inc ~~ RI_gpd
  RI_inc ~~ 0*inc_w_t1
  RI_inc ~~ 0*gpd_w_t1
  RI_gpd ~~ 0*inc_w_t1
  RI_gpd ~~ 0*gpd_w_t1
'

########## Standardized effects with 95% CIs ###################################

std_ci_table <- function(fit, intercepts = c("RI_inc", "RI_gpd")) {
  std_results <- standardizedSolution(fit, type = "std.all",
                                      ci = TRUE, level = 0.95)
  
  keep <- std_results$op == "~" |
    (std_results$op == "~~" &
       std_results$lhs == intercepts[1] &
       std_results$rhs == intercepts[2])
  
  out <- std_results[keep,
                     c("lhs", "op", "rhs", "est.std", "se", "ci.lower", "ci.upper")]
  rownames(out) <- NULL
  out
}

print_std_ci <- function(fit, intercepts = c("RI_inc", "RI_gpd")) {
  cat("\n\nStandardized (std.all) regression paths and random-intercept\n")
  cat("correlation (", intercepts[1], " ~~ ", intercepts[2], "), with 95% ",
      "confidence intervals:\n\n", sep = "")
  if (!lavInspect(fit, "converged")) {
    cat("WARNING: this model did not converge, so the values below are ",
        "unreliable (NA = not computable).\n\n", sep = "")
  }
  print(std_ci_table(fit, intercepts))
}

########## Bern: unconstrained model #############################################

fit_riclpm_bern <- sem(
  model_riclpm_4w,
  data = df_bern,
  missing = "ML")

summary(fit_riclpm_bern, standardized = TRUE, fit.measures = TRUE, ci = TRUE)

std_ci_table(fit_riclpm_bern)

########## Bern: constrained model ################################################

fit_riclpm_bern_constrained <- sem(
  model_riclpm_4w_constrained,
  data = df_bern,
  missing = "ML")

summary(fit_riclpm_bern_constrained, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_riclpm_bern_constrained)

########## Bern: constrained vs. unconstrained ####################################

summary(compareFit(fit_riclpm_bern_constrained, fit_riclpm_bern))

fitMeasures(
  fit_riclpm_bern_constrained,
  "rmsea",
  h1.model = fit_riclpm_bern)

########## Zurich: unconstrained model ############################################

fit_riclpm_zurich <- sem(
  model_riclpm_4w,
  data = df_zurich,
  missing = "ML")

summary(fit_riclpm_zurich, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_riclpm_zurich)

########## Zurich: constrained model ##############################################

fit_riclpm_zurich_constrained <- sem(
  model_riclpm_4w_constrained,
  data = df_zurich,
  missing = "ML")

summary(fit_riclpm_zurich_constrained, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_riclpm_zurich_constrained)

########## Zurich: constrained vs. unconstrained ##################################

summary(compareFit(fit_riclpm_zurich_constrained, fit_riclpm_zurich))

fitMeasures(
  fit_riclpm_zurich_constrained,
  "rmsea",
  h1.model = fit_riclpm_zurich)

########## Trier: unconstrained model, 6 waves ####################################

fit_riclpm_trier <- sem(
  model_riclpm_6w,
  data = df_trier,
  missing = "ML")

summary(fit_riclpm_trier, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_riclpm_trier)

########## Trier: constrained model, 6 waves ######################################

fit_riclpm_trier_constrained <- sem(
  model_riclpm_6w_constrained,
  data = df_trier,
  missing = "ML")

summary(fit_riclpm_trier_constrained, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_riclpm_trier_constrained)

########## Trier: constrained vs. unconstrained ###################################

summary(compareFit(fit_riclpm_trier_constrained, fit_riclpm_trier))

fitMeasures(
  fit_riclpm_trier_constrained,
  "rmsea",
  h1.model = fit_riclpm_trier)

########## Sample size and convergence check ######################################

riclpm_sample_check <- tibble(
  sample = c("Bern", "Zurich", "Trier"),
  unconstrained_n = c(
    lavInspect(fit_riclpm_bern, "nobs"),
    lavInspect(fit_riclpm_zurich, "nobs"),
    lavInspect(fit_riclpm_trier, "nobs")),
  constrained_n = c(
    lavInspect(fit_riclpm_bern_constrained, "nobs"),
    lavInspect(fit_riclpm_zurich_constrained, "nobs"),
    lavInspect(fit_riclpm_trier_constrained, "nobs")),
  expected_n = c(101L, 145L, 1410L),
  unconstrained_converged = c(
    lavInspect(fit_riclpm_bern, "converged"),
    lavInspect(fit_riclpm_zurich, "converged"),
    lavInspect(fit_riclpm_trier, "converged")),
  constrained_converged = c(
    lavInspect(fit_riclpm_bern_constrained, "converged"),
    lavInspect(fit_riclpm_zurich_constrained, "converged"),
    lavInspect(fit_riclpm_trier_constrained, "converged"))
)

print(riclpm_sample_check)

if (
  any(riclpm_sample_check$unconstrained_n != riclpm_sample_check$expected_n) ||
  any(riclpm_sample_check$constrained_n != riclpm_sample_check$expected_n) ||
  any(!riclpm_sample_check$unconstrained_converged) ||
  any(!riclpm_sample_check$constrained_converged)
) {
  stop("At least one RI-CLPM failed its sample-size or convergence check.")
}

################################################################################
########## Save model output as text files ####################################
#
# Each model's raw R/lavaan output (the unconstrained summary, the
# constrained summary, and the constrained-vs-unconstrained comparison),
# exactly as summary()/compareFit()/fitMeasures() produce it above. Each
# summary file ends with the standardized paths' 95% CI table.

dir.create("results/03 riclpm/model summaries", recursive = TRUE, showWarnings = FALSE)
dir.create("results/03 riclpm/model comparisons", recursive = TRUE, showWarnings = FALSE)

save_output <- function(file, expr) {
  writeLines(capture.output(force(expr)), con = file)
}

save_output(
  "results/03 riclpm/model summaries/bern_4w_unconstrained_summary.txt",
  {
    print(summary(fit_riclpm_bern, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_riclpm_bern)
  }
)

save_output(
  "results/03 riclpm/model summaries/bern_4w_constrained_summary.txt",
  {
    print(summary(fit_riclpm_bern_constrained, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_riclpm_bern_constrained)
  }
)

save_output(
  "results/03 riclpm/model comparisons/bern_4w_model_comparison.txt",
  {
    print(summary(compareFit(fit_riclpm_bern_constrained, fit_riclpm_bern)))
    cat("\n\nRMSEA using unconstrained model as h1 model:\n")
    print(fitMeasures(
      fit_riclpm_bern_constrained,
      "rmsea",
      h1.model = fit_riclpm_bern
    ))
  }
)

save_output(
  "results/03 riclpm/model summaries/zurich_4w_unconstrained_summary.txt",
  {
    print(summary(fit_riclpm_zurich, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_riclpm_zurich)
  }
)

save_output(
  "results/03 riclpm/model summaries/zurich_4w_constrained_summary.txt",
  {
    print(summary(fit_riclpm_zurich_constrained, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_riclpm_zurich_constrained)
  }
)

save_output(
  "results/03 riclpm/model comparisons/zurich_4w_model_comparison.txt",
  {
    print(summary(compareFit(fit_riclpm_zurich_constrained, fit_riclpm_zurich)))
    cat("\n\nRMSEA using unconstrained model as h1 model:\n")
    print(fitMeasures(
      fit_riclpm_zurich_constrained,
      "rmsea",
      h1.model = fit_riclpm_zurich
    ))
  }
)

save_output(
  "results/03 riclpm/model summaries/trier_6w_unconstrained_summary.txt",
  {
    print(summary(fit_riclpm_trier, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_riclpm_trier)
  }
)

save_output(
  "results/03 riclpm/model summaries/trier_6w_constrained_summary.txt",
  {
    print(summary(fit_riclpm_trier_constrained, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_riclpm_trier_constrained)
  }
)

save_output(
  "results/03 riclpm/model comparisons/trier_6w_model_comparison.txt",
  {
    print(summary(compareFit(fit_riclpm_trier_constrained, fit_riclpm_trier)))
    cat("\n\nRMSEA using unconstrained model as h1 model:\n")
    print(fitMeasures(
      fit_riclpm_trier_constrained,
      "rmsea",
      h1.model = fit_riclpm_trier
    ))
  }
)

################################################################################
########## Plot the constrained RI-CLPMs ######################################

# Layout for the 4-wave models (Bern and Zurich): random intercepts on top
# and bottom, observed scores next to the random intercepts, wave-specific
# residual factors in the middle rows.
lay_riclpm_4w <- get_layout(
  NA,         NA,         "RI_inc",   NA,
  "ink1",     "ink2",     "ink3",     "ink4",
  "inc_w_t1", "inc_w_t2", "inc_w_t3", "inc_w_t4",
  "gpd_w_t1", "gpd_w_t2", "gpd_w_t3", "gpd_w_t4",
  "gpd1",     "gpd2",     "gpd3",     "gpd4",
  NA,         NA,         "RI_gpd",   NA,
  rows = 6)

# Same layout, extended to six waves for Trier.
lay_riclpm_6w <- get_layout(
  NA,         NA,         "RI_inc",   NA,         NA,         NA,
  "ink1",     "ink2",     "ink3",     "ink4",     "ink5",     "ink6",
  "inc_w_t1", "inc_w_t2", "inc_w_t3", "inc_w_t4", "inc_w_t5", "inc_w_t6",
  "gpd_w_t1", "gpd_w_t2", "gpd_w_t3", "gpd_w_t4", "gpd_w_t5", "gpd_w_t6",
  "gpd1",     "gpd2",     "gpd3",     "gpd4",     "gpd5",     "gpd6",
  NA,         NA,         "RI_gpd",   NA,         NA,         NA,
  rows = 6)

########## Bern plots #############################################################

p_riclpm_bern_names <- graph_sem(
  fit_riclpm_bern_constrained,
  layout = lay_riclpm_4w,
  label = "name")

ggsave(
  "plots/03 riclpm/riclpm_bern_names.png",
  plot = p_riclpm_bern_names,
  width = 12, height = 10, dpi = 300, bg = "white")

p_riclpm_bern_est <- graph_sem(
  fit_riclpm_bern_constrained,
  layout = lay_riclpm_4w,
  label = "est_std")

ggsave(
  "plots/03 riclpm/riclpm_bern_est_std.png",
  plot = p_riclpm_bern_est,
  width = 12, height = 10, dpi = 300, bg = "white")

########## Zurich plots ############################################################

p_riclpm_zurich_names <- graph_sem(
  fit_riclpm_zurich_constrained,
  layout = lay_riclpm_4w,
  label = "name")

ggsave(
  "plots/03 riclpm/riclpm_zurich_names.png",
  plot = p_riclpm_zurich_names,
  width = 12, height = 10, dpi = 300, bg = "white")

p_riclpm_zurich_est <- graph_sem(
  fit_riclpm_zurich_constrained,
  layout = lay_riclpm_4w,
  label = "est_std")

ggsave(
  "plots/03 riclpm/riclpm_zurich_est_std.png",
  plot = p_riclpm_zurich_est,
  width = 12, height = 10, dpi = 300, bg = "white")

########## Trier plots #############################################################

p_riclpm_trier_names <- graph_sem(
  fit_riclpm_trier_constrained,
  layout = lay_riclpm_6w,
  label = "name")

ggsave(
  "plots/03 riclpm/riclpm_trier_names.png",
  plot = p_riclpm_trier_names,
  width = 14, height = 10, dpi = 300, bg = "white")

p_riclpm_trier_est <- graph_sem(
  fit_riclpm_trier_constrained,
  layout = lay_riclpm_6w,
  label = "est_std")

ggsave(
  "plots/03 riclpm/riclpm_trier_est_std.png",
  plot = p_riclpm_trier_est,
  width = 14, height = 10, dpi = 300, bg = "white")
