# Step 1 of 4: tests whether the raw wave means show a structural change
# (a trend) over time. 

library(readr)
library(dplyr)
library(psych)
library(tibble)
library(stringr)
library(lavaan)
library(semTools)
library(ggplot2)

dir.create("plots/01 mean structure", recursive = TRUE, showWarnings = FALSE)

########## Data ####################################################################

load("data/data.RData")

########## Model syntax ############################################################

# All variances and covariances among a set of observed variables, written
# out explicitly so the covariance structure is genuinely saturated (lavaan
# does not add these automatically for a set of plain observed variables
# that only have an intercept formula).
saturated_cov_syntax <- function(vars) {
  variances <- paste0(vars, " ~~ ", vars)
  covariances <- combn(vars, 2, function(pair) paste0(pair[1], " ~~ ", pair[2]))
  paste(c(variances, covariances), collapse = "\n")
}

vars_4w <- c("ink1", "ink2", "ink3", "ink4", "gpd1", "gpd2", "gpd3", "gpd4")
vars_6w <- c("ink1", "ink2", "ink3", "ink4", "ink5", "ink6", "gpd1", "gpd2", "gpd3", "gpd4", "gpd5", "gpd6")

# Means constrained equal across waves, one shared mean per construct.
model_meanstructure_4w <- paste0('
  ink1 ~ mu_inc*1
  ink2 ~ mu_inc*1
  ink3 ~ mu_inc*1
  ink4 ~ mu_inc*1

  gpd1 ~ mu_gpd*1
  gpd2 ~ mu_gpd*1
  gpd3 ~ mu_gpd*1
  gpd4 ~ mu_gpd*1
', saturated_cov_syntax(vars_4w))

model_meanstructure_6w <- paste0('
  ink1 ~ mu_inc*1
  ink2 ~ mu_inc*1
  ink3 ~ mu_inc*1
  ink4 ~ mu_inc*1
  ink5 ~ mu_inc*1
  ink6 ~ mu_inc*1

  gpd1 ~ mu_gpd*1
  gpd2 ~ mu_gpd*1
  gpd3 ~ mu_gpd*1
  gpd4 ~ mu_gpd*1
  gpd5 ~ mu_gpd*1
  gpd6 ~ mu_gpd*1
', saturated_cov_syntax(vars_6w))

# Companion model: every wave's mean estimated freely (same saturated
# covariance structure), used only to plot the observed trend.
model_freemeans_4w <- paste0('
  ink1 ~ 1
  ink2 ~ 1
  ink3 ~ 1
  ink4 ~ 1

  gpd1 ~ 1
  gpd2 ~ 1
  gpd3 ~ 1
  gpd4 ~ 1
', saturated_cov_syntax(vars_4w))

model_freemeans_6w <- paste0('
  ink1 ~ 1
  ink2 ~ 1
  ink3 ~ 1
  ink4 ~ 1
  ink5 ~ 1
  ink6 ~ 1

  gpd1 ~ 1
  gpd2 ~ 1
  gpd3 ~ 1
  gpd4 ~ 1
  gpd5 ~ 1
  gpd6 ~ 1
', saturated_cov_syntax(vars_6w))

########## Bern ####################################################################

fit_meanstructure_bern <- sem(
  model_meanstructure_4w,
  data = df_bern,
  missing = "ML")

summary(fit_meanstructure_bern, fit.measures = TRUE)

fit_freemeans_bern <- sem(
  model_freemeans_4w,
  data = df_bern,
  missing = "ML")

fitMeasures(
  fit_meanstructure_bern,
  "rmsea",
  h1.model = fit_freemeans_bern)

########## Zurich ##################################################################

fit_meanstructure_zurich <- sem(
  model_meanstructure_4w,
  data = df_zurich,
  missing = "ML")

summary(fit_meanstructure_zurich, fit.measures = TRUE)

fit_freemeans_zurich <- sem(
  model_freemeans_4w,
  data = df_zurich,
  missing = "ML")

fitMeasures(
  fit_meanstructure_zurich,
  "rmsea",
  h1.model = fit_freemeans_zurich)

########## Trier, 6 waves ##########################################################

fit_meanstructure_trier <- sem(
  model_meanstructure_6w,
  data = df_trier,
  missing = "ML")

summary(fit_meanstructure_trier, fit.measures = TRUE)

fit_freemeans_trier <- sem(
  model_freemeans_6w,
  data = df_trier,
  missing = "ML")

fitMeasures(
  fit_meanstructure_trier,
  "rmsea",
  h1.model = fit_freemeans_trier)

########## Sample size and convergence check #######################################

meanstructure_sample_check <- tibble(
  sample = c("Bern", "Zurich", "Trier"),
  n = c(
    lavInspect(fit_meanstructure_bern, "nobs"),
    lavInspect(fit_meanstructure_zurich, "nobs"),
    lavInspect(fit_meanstructure_trier, "nobs")),
  expected_n = c(101L, 145L, 1410L),
  meanstructure_converged = c(
    lavInspect(fit_meanstructure_bern, "converged"),
    lavInspect(fit_meanstructure_zurich, "converged"),
    lavInspect(fit_meanstructure_trier, "converged")),
  freemeans_converged = c(
    lavInspect(fit_freemeans_bern, "converged"),
    lavInspect(fit_freemeans_zurich, "converged"),
    lavInspect(fit_freemeans_trier, "converged"))
)

print(meanstructure_sample_check)

if (
  any(meanstructure_sample_check$n != meanstructure_sample_check$expected_n) ||
  any(!meanstructure_sample_check$meanstructure_converged) ||
  any(!meanstructure_sample_check$freemeans_converged)
) {
  stop("At least one mean-structure model failed its sample-size or convergence check.")
}

################################################################################
########## Save model output as text files #####################################
#
# Each sample's raw R/lavaan output: the mean-constrained model's summary
# (its chi-square test is the structural-change test itself) and a plain
# verdict line.

dir.create("results/01 mean structure", recursive = TRUE, showWarnings = FALSE)

save_output <- function(file, expr) {
  writeLines(capture.output(force(expr)), con = file)
}

save_meanstructure_result <- function(file, fit, sample_name) {
  save_output(file, {
    print(summary(fit, fit.measures = TRUE))
    cat("\n\nStructural-change test (Hamaker et al., 2015): means constrained\n")
    cat("equal across waves, covariance structure saturated. A significant\n")
    cat("chi-square means the wave means genuinely differ over time.\n\n")
    chisq <- fitMeasures(fit, "chisq")
    df <- fitMeasures(fit, "df")
    pvalue <- fitMeasures(fit, "pvalue")
    verdict <- if (pvalue < .05) {
      "there IS a significant structural change (trend) over time."
    } else {
      "there is NO significant structural change (trend) over time."
    }
    cat(sample_name, ": chi-square(", df, ") = ", round(chisq, 2),
        ", p = ", format.pval(pvalue, digits = 3), " -- ", verdict, "\n",
        sep = "")
  })
}

save_meanstructure_result(
  "results/01 mean structure/bern_4w_mean_structure.txt",
  fit_meanstructure_bern, "Bern")

save_meanstructure_result(
  "results/01 mean structure/zurich_4w_mean_structure.txt",
  fit_meanstructure_zurich, "Zurich")

save_meanstructure_result(
  "results/01 mean structure/trier_6w_mean_structure.txt",
  fit_meanstructure_trier, "Trier")

################################################################################
########## Plot the wave means #########################################

# Uses the freely estimated means (same saturated covariance structure) to
# show each construct's trajectory, with the FIML-based standard
# errors as 95% confidence bars.
wave_means <- function(fit, sample_name, waves) {
  parameterEstimates(fit) %>%
    filter(op == "~1") %>%
    mutate(
      sample = sample_name,
      construct = ifelse(str_starts(lhs, "ink"), "INC", "GPD"),
      wave = rep(seq_len(waves), times = 2)
    ) %>%
    select(sample, construct, wave, est, se) %>%
    as_tibble()
}

mean_trends <- bind_rows(
  wave_means(fit_freemeans_bern, "Bern", 4),
  wave_means(fit_freemeans_zurich, "Zurich", 4),
  wave_means(fit_freemeans_trier, "Trier", 6)
)

print(mean_trends, n = Inf)

plot_mean_trend <- function(data, sample_name) {
  ggplot(filter(data, sample == sample_name), aes(x = wave, y = est)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), width = 0.15) +
    facet_wrap(~construct, scales = "free_y") +
    scale_x_continuous(breaks = seq_len(max(data$wave[data$sample == sample_name]))) +
    labs(
      title = paste0(sample_name, ": wave means (95% CI)"),
      x = "Wave", y = "Mean") +
    theme_minimal(base_size = 13)
}

ggsave(
  "plots/01 mean structure/mean_trend_bern.png",
  plot = plot_mean_trend(mean_trends, "Bern"),
  width = 9, height = 5, dpi = 300, bg = "white")

ggsave(
  "plots/01 mean structure/mean_trend_zurich.png",
  plot = plot_mean_trend(mean_trends, "Zurich"),
  width = 9, height = 5, dpi = 300, bg = "white")

ggsave(
  "plots/01 mean structure/mean_trend_trier.png",
  plot = plot_mean_trend(mean_trends, "Trier"),
  width = 9, height = 5, dpi = 300, bg = "white")
