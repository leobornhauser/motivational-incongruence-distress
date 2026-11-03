# Step 2 of 4: compares a nondetrended RI-CLPM (intercepts constrained
# equal across waves) against the detrended RI-CLPM (Falkenstrom et al.,
# 2023) already fitted in 03 riclpm/, which this script sources directly

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

dir.create("plots/02 riclpm mean structure", recursive = TRUE, showWarnings = FALSE)

########## The detrended RI-CLPM (step 3) #########################################

# Provides df_bern/df_zurich/df_trier, the fitted detrended models, and
# lay_riclpm_4w/6w for plotting.
source("scripts/models/03_riclpm.R")

########## Define the nondetrended RI-CLPM: 4 waves (Bern and Zurich) ############

model_riclpm_nondetrended_4w <- '

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

  # nondetrended: one shared intercept per construct, held equal over waves
  ink1 ~ c_inc*1
  ink2 ~ c_inc*1
  ink3 ~ c_inc*1
  ink4 ~ c_inc*1

  gpd1 ~ c_gpd*1
  gpd2 ~ c_gpd*1
  gpd3 ~ c_gpd*1
  gpd4 ~ c_gpd*1

  # structural relations between residual factors
  inc_w_t2 ~ inc_w_t1
  inc_w_t3 ~ inc_w_t2
  inc_w_t4 ~ inc_w_t3

  gpd_w_t2 ~ gpd_w_t1
  gpd_w_t3 ~ gpd_w_t2
  gpd_w_t4 ~ gpd_w_t3

  gpd_w_t2 ~ inc_w_t1
  gpd_w_t3 ~ inc_w_t2
  gpd_w_t4 ~ inc_w_t3

  inc_w_t2 ~ gpd_w_t1
  inc_w_t3 ~ gpd_w_t2
  inc_w_t4 ~ gpd_w_t3

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

model_riclpm_nondetrended_4w_constrained <- '

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

  ink1 ~ c_inc*1
  ink2 ~ c_inc*1
  ink3 ~ c_inc*1
  ink4 ~ c_inc*1

  gpd1 ~ c_gpd*1
  gpd2 ~ c_gpd*1
  gpd3 ~ c_gpd*1
  gpd4 ~ c_gpd*1

  # structural relations between residual factors
  inc_w_t2 ~ a*inc_w_t1
  inc_w_t3 ~ a*inc_w_t2
  inc_w_t4 ~ a*inc_w_t3

  gpd_w_t2 ~ b*gpd_w_t1
  gpd_w_t3 ~ b*gpd_w_t2
  gpd_w_t4 ~ b*gpd_w_t3

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

########## Define the nondetrended RI-CLPM: 6 waves (Trier) #######################

model_riclpm_nondetrended_6w <- '

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

  ink1 ~ c_inc*1
  ink2 ~ c_inc*1
  ink3 ~ c_inc*1
  ink4 ~ c_inc*1
  ink5 ~ c_inc*1
  ink6 ~ c_inc*1

  gpd1 ~ c_gpd*1
  gpd2 ~ c_gpd*1
  gpd3 ~ c_gpd*1
  gpd4 ~ c_gpd*1
  gpd5 ~ c_gpd*1
  gpd6 ~ c_gpd*1

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

model_riclpm_nondetrended_6w_constrained <- '

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

  ink1 ~ c_inc*1
  ink2 ~ c_inc*1
  ink3 ~ c_inc*1
  ink4 ~ c_inc*1
  ink5 ~ c_inc*1
  ink6 ~ c_inc*1

  gpd1 ~ c_gpd*1
  gpd2 ~ c_gpd*1
  gpd3 ~ c_gpd*1
  gpd4 ~ c_gpd*1
  gpd5 ~ c_gpd*1
  gpd6 ~ c_gpd*1

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

# The nondetrended, unconstrained combination sometimes fails to converge
# (see Bern below); comparisons involving it are guarded rather than
# stopping the whole script.

print_comparison <- function(fit_constrained, fit_unconstrained, sample_name) {
  if (lavInspect(fit_unconstrained, "converged")) {
    print(summary(compareFit(fit_constrained, fit_unconstrained)))
    cat("\n\nRMSEA using unconstrained model as h1 model:\n")
    print(fitMeasures(fit_constrained, "rmsea", h1.model = fit_unconstrained))
  } else {
    message(
      sample_name, ": the unconstrained nondetrended model did not ",
      "converge, so the constrained vs. unconstrained comparison is skipped.")
  }
}

# The nondetrended model is nested within the detrended model, so this is
# a nested chi-square test of the detrending question itself.
print_detrending_comparison <- function(fit_nondetrended_constrained, fit_detrended_constrained, sample_name) {
  if (lavInspect(fit_nondetrended_constrained, "converged")) {
    print(summary(compareFit(fit_nondetrended_constrained, fit_detrended_constrained)))
    cat("\n\nRMSEA using detrended model as h1 model:\n")
    print(fitMeasures(fit_nondetrended_constrained, "rmsea", h1.model = fit_detrended_constrained))
  } else {
    message(sample_name, ": the detrending comparison is skipped because the nondetrended constrained model did not converge.")
  }
}

########## Bern #####################################################################

fit_riclpm_nondetrended_bern <- sem(
  model_riclpm_nondetrended_4w,
  data = df_bern,
  missing = "ML")

summary(fit_riclpm_nondetrended_bern, standardized = TRUE, fit.measures = TRUE)

fit_riclpm_nondetrended_bern_constrained <- sem(
  model_riclpm_nondetrended_4w_constrained,
  data = df_bern,
  missing = "ML")

summary(fit_riclpm_nondetrended_bern_constrained, standardized = TRUE, fit.measures = TRUE)

print_comparison(fit_riclpm_nondetrended_bern_constrained, fit_riclpm_nondetrended_bern, "Bern")
print_detrending_comparison(fit_riclpm_nondetrended_bern_constrained, fit_riclpm_bern_constrained, "Bern")

########## Zurich ####################################################################

fit_riclpm_nondetrended_zurich <- sem(
  model_riclpm_nondetrended_4w,
  data = df_zurich,
  missing = "ML")

summary(fit_riclpm_nondetrended_zurich, standardized = TRUE, fit.measures = TRUE)

fit_riclpm_nondetrended_zurich_constrained <- sem(
  model_riclpm_nondetrended_4w_constrained,
  data = df_zurich,
  missing = "ML")

summary(fit_riclpm_nondetrended_zurich_constrained, standardized = TRUE, fit.measures = TRUE)

print_comparison(fit_riclpm_nondetrended_zurich_constrained, fit_riclpm_nondetrended_zurich, "Zurich")
print_detrending_comparison(fit_riclpm_nondetrended_zurich_constrained, fit_riclpm_zurich_constrained, "Zurich")

########## Trier, 6 waves ############################################################

fit_riclpm_nondetrended_trier <- sem(
  model_riclpm_nondetrended_6w,
  data = df_trier,
  missing = "ML")

summary(fit_riclpm_nondetrended_trier, standardized = TRUE, fit.measures = TRUE)

fit_riclpm_nondetrended_trier_constrained <- sem(
  model_riclpm_nondetrended_6w_constrained,
  data = df_trier,
  missing = "ML")

summary(fit_riclpm_nondetrended_trier_constrained, standardized = TRUE, fit.measures = TRUE)

print_comparison(fit_riclpm_nondetrended_trier_constrained, fit_riclpm_nondetrended_trier, "Trier")
print_detrending_comparison(fit_riclpm_nondetrended_trier_constrained, fit_riclpm_trier_constrained, "Trier")

########## Sample size and convergence check #########################################

# The detrended models were already checked by 03 riclpm/'s own check when
# it was sourced above; only the new nondetrended models need checking here.
riclpm_nondetrended_sample_check <- tibble(
  sample = rep(c("Bern", "Zurich", "Trier"), each = 2),
  model = rep(c("unconstrained", "constrained"), times = 3),
  n = c(
    lavInspect(fit_riclpm_nondetrended_bern, "nobs"),
    lavInspect(fit_riclpm_nondetrended_bern_constrained, "nobs"),
    lavInspect(fit_riclpm_nondetrended_zurich, "nobs"),
    lavInspect(fit_riclpm_nondetrended_zurich_constrained, "nobs"),
    lavInspect(fit_riclpm_nondetrended_trier, "nobs"),
    lavInspect(fit_riclpm_nondetrended_trier_constrained, "nobs")),
  expected_n = rep(c(101L, 145L, 1410L), each = 2),
  converged = c(
    lavInspect(fit_riclpm_nondetrended_bern, "converged"),
    lavInspect(fit_riclpm_nondetrended_bern_constrained, "converged"),
    lavInspect(fit_riclpm_nondetrended_zurich, "converged"),
    lavInspect(fit_riclpm_nondetrended_zurich_constrained, "converged"),
    lavInspect(fit_riclpm_nondetrended_trier, "converged"),
    lavInspect(fit_riclpm_nondetrended_trier_constrained, "converged"))
)

print(riclpm_nondetrended_sample_check, n = Inf)

# Only the constrained model needs to converge (it is what gets plotted);
# an unconstrained failure (Bern) is already reported via message() above.
riclpm_nondetrended_constrained_check <- riclpm_nondetrended_sample_check %>%
  filter(model == "constrained")

if (
  any(riclpm_nondetrended_constrained_check$n != riclpm_nondetrended_constrained_check$expected_n) ||
  any(!riclpm_nondetrended_constrained_check$converged)
) {
  stop("At least one constrained nondetrended RI-CLPM failed its sample-size or convergence check.")
}

if (any(!riclpm_nondetrended_sample_check$converged)) {
  message(
    "Note: the unconstrained nondetrended model did not converge for at ",
    "least one sample -- see the messages above, next to each sample's ",
    "within-model comparison. The constrained model is unaffected and is ",
    "what is plotted and used for the detrending question above.")
}

################################################################################
########## Save model output as text files ####################################
#
# Only the nondetrended models: the detrended models' output already lives
# in results/03 riclpm/.

summaries_dir <- "results/02 riclpm mean structure/model summaries"
comparisons_dir <- "results/02 riclpm mean structure/model comparisons"
dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparisons_dir, recursive = TRUE, showWarnings = FALSE)

save_output <- function(file, expr) {
  writeLines(capture.output(force(expr)), con = file)
}

# Mirrors print_comparison() above, writing to a file instead of the console.
save_model_comparison <- function(file, fit_constrained, fit_unconstrained, sample_name) {
  if (lavInspect(fit_unconstrained, "converged")) {
    save_output(file, {
      print(summary(compareFit(fit_constrained, fit_unconstrained)))
      cat("\n\nRMSEA using unconstrained model as h1 model:\n")
      print(fitMeasures(fit_constrained, "rmsea", h1.model = fit_unconstrained))
    })
  } else {
    save_output(file, {
      cat(
        sample_name, ": the unconstrained nondetrended model did not ",
        "converge, so the constrained vs. unconstrained comparison could ",
        "not be computed. See the constrained model's own summary file ",
        "instead.\n",
        sep = ""
      )
    })
  }
}

save_detrending_comparison <- function(file, fit_nondetrended_constrained, fit_detrended_constrained, sample_name) {
  if (lavInspect(fit_nondetrended_constrained, "converged")) {
    save_output(file, {
      print(summary(compareFit(fit_nondetrended_constrained, fit_detrended_constrained)))
      cat("\n\nRMSEA using detrended model as h1 model:\n")
      print(fitMeasures(fit_nondetrended_constrained, "rmsea", h1.model = fit_detrended_constrained))
    })
  } else {
    save_output(file, {
      cat(
        sample_name, ": the detrending comparison could not be computed ",
        "because the nondetrended constrained model did not converge.\n",
        sep = ""
      )
    })
  }
}

for (sample_name in c("bern", "zurich", "trier")) {
  waves <- if (sample_name == "trier") "6w" else "4w"
  sample_label <- str_to_title(sample_name)
  fit_nd <- get(paste0("fit_riclpm_nondetrended_", sample_name))
  fit_nd_c <- get(paste0("fit_riclpm_nondetrended_", sample_name, "_constrained"))
  fit_d_c <- get(paste0("fit_riclpm_", sample_name, "_constrained"))

  save_output(
    file.path(summaries_dir, paste0("nondetrended_", sample_name, "_", waves, "_unconstrained_summary.txt")),
    summary(fit_nd, standardized = TRUE, fit.measures = TRUE)
  )
  save_output(
    file.path(summaries_dir, paste0("nondetrended_", sample_name, "_", waves, "_constrained_summary.txt")),
    summary(fit_nd_c, standardized = TRUE, fit.measures = TRUE)
  )
  save_model_comparison(
    file.path(comparisons_dir, paste0("nondetrended_", sample_name, "_", waves, "_model_comparison.txt")),
    fit_nd_c, fit_nd, sample_label
  )
  save_detrending_comparison(
    file.path(comparisons_dir, paste0(sample_name, "_", waves, "_detrending_comparison.txt")),
    fit_nd_c, fit_d_c, sample_label
  )
}

################################################################################
########## Plot the constrained nondetrended RI-CLPMs #########################
#
# The detrended path diagrams are already in plots/03 riclpm/; lay_riclpm_4w/6w
# come from the source() call above.

########## Bern plot #################################################################

p_nondetrended_bern_names <- graph_sem(
  fit_riclpm_nondetrended_bern_constrained,
  layout = lay_riclpm_4w,
  label = "name")

ggsave(
  "plots/02 riclpm mean structure/riclpm_nondetrended_bern_names.png",
  plot = p_nondetrended_bern_names,
  width = 12, height = 10, dpi = 300, bg = "white")

p_nondetrended_bern_est <- graph_sem(
  fit_riclpm_nondetrended_bern_constrained,
  layout = lay_riclpm_4w,
  label = "est_std")

ggsave(
  "plots/02 riclpm mean structure/riclpm_nondetrended_bern_est_std.png",
  plot = p_nondetrended_bern_est,
  width = 12, height = 10, dpi = 300, bg = "white")

########## Zurich plot ################################################################

p_nondetrended_zurich_names <- graph_sem(
  fit_riclpm_nondetrended_zurich_constrained,
  layout = lay_riclpm_4w,
  label = "name")

ggsave(
  "plots/02 riclpm mean structure/riclpm_nondetrended_zurich_names.png",
  plot = p_nondetrended_zurich_names,
  width = 12, height = 10, dpi = 300, bg = "white")

p_nondetrended_zurich_est <- graph_sem(
  fit_riclpm_nondetrended_zurich_constrained,
  layout = lay_riclpm_4w,
  label = "est_std")

ggsave(
  "plots/02 riclpm mean structure/riclpm_nondetrended_zurich_est_std.png",
  plot = p_nondetrended_zurich_est,
  width = 12, height = 10, dpi = 300, bg = "white")

########## Trier plot #################################################################

p_nondetrended_trier_names <- graph_sem(
  fit_riclpm_nondetrended_trier_constrained,
  layout = lay_riclpm_6w,
  label = "name")

ggsave(
  "plots/02 riclpm mean structure/riclpm_nondetrended_trier_names.png",
  plot = p_nondetrended_trier_names,
  width = 14, height = 10, dpi = 300, bg = "white")

p_nondetrended_trier_est <- graph_sem(
  fit_riclpm_nondetrended_trier_constrained,
  layout = lay_riclpm_6w,
  label = "est_std")

ggsave(
  "plots/02 riclpm mean structure/riclpm_nondetrended_trier_est_std.png",
  plot = p_nondetrended_trier_est,
  width = 14, height = 10, dpi = 300, bg = "white")
