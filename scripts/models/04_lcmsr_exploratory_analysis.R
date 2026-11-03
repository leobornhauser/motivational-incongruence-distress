# Step 4 of 4: exploratory check on step 3, using a bivariate latent curve
# model with structured residuals (LCM-SR; Curran et al., 2014).

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

dir.create(
  "plots/04 lcm-sr exploratory analysis",
  recursive = TRUE,
  showWarnings = FALSE
)
########## Data ###################################################################

load("data/data.RData")

########## Define the LCM-SR: 4 waves (Bern and Zurich) #########################

# Unconstrained: autoregressive and cross-lagged paths free to differ across
# waves. The slope loadings 0, 1, 2, 3 index the four nominal assessment
# waves (actual assessment timing varied somewhat within each sample, so a
# slope is best read as average change per assessment interval).
model_lcmsr_4w <- '

  # random-intercept and random-slope factors: a random intercept (loading
  # of 1 on every wave, exactly as in the main RI-CLPM) and a random linear
  # slope (loadings 0, 1, 2, 3)
  I_inc =~ 1*ink1 + 1*ink2 + 1*ink3 + 1*ink4
  S_inc =~ 0*ink1 + 1*ink2 + 2*ink3 + 3*ink4

  I_gpd =~ 1*gpd1 + 1*gpd2 + 1*gpd3 + 1*gpd4
  S_gpd =~ 0*gpd1 + 1*gpd2 + 2*gpd3 + 3*gpd4

  # wave-specific residual factors
  inc_w_t1 =~ 1*ink1
  inc_w_t2 =~ 1*ink2
  inc_w_t3 =~ 1*ink3
  inc_w_t4 =~ 1*ink4

  gpd_w_t1 =~ 1*gpd1
  gpd_w_t2 =~ 1*gpd2
  gpd_w_t3 =~ 1*gpd3
  gpd_w_t4 =~ 1*gpd4

  # the observed scores carry no residual variance or mean of their own
  ink1 ~~ 0*ink1
  ink2 ~~ 0*ink2
  ink3 ~~ 0*ink3
  ink4 ~~ 0*ink4
  ink1 ~ 0*1
  ink2 ~ 0*1
  ink3 ~ 0*1
  ink4 ~ 0*1

  gpd1 ~~ 0*gpd1
  gpd2 ~~ 0*gpd2
  gpd3 ~~ 0*gpd3
  gpd4 ~~ 0*gpd4
  gpd1 ~ 0*1
  gpd2 ~ 0*1
  gpd3 ~ 0*1
  gpd4 ~ 0*1

  # the wave-specific residual factors carry no mean of their own
  # all systematic change over time is carried by the growth factor means
  inc_w_t1 ~ 0*1
  inc_w_t2 ~ 0*1
  inc_w_t3 ~ 0*1
  inc_w_t4 ~ 0*1

  gpd_w_t1 ~ 0*1
  gpd_w_t2 ~ 0*1
  gpd_w_t3 ~ 0*1
  gpd_w_t4 ~ 0*1

  # growth factor means
  I_inc ~ 1
  S_inc ~ 1
  I_gpd ~ 1
  S_gpd ~ 1

  # the two intercepts and two slopes are all free to vary and to covary
  # with each other
  I_inc ~~ I_inc
  S_inc ~~ S_inc
  I_gpd ~~ I_gpd
  S_gpd ~~ S_gpd

  I_inc ~~ S_inc
  I_inc ~~ I_gpd
  I_inc ~~ S_gpd
  S_inc ~~ I_gpd
  S_inc ~~ S_gpd
  I_gpd ~~ S_gpd

  # wave-specific state variances
  inc_w_t1 ~~ inc_w_t1
  inc_w_t2 ~~ inc_w_t2
  inc_w_t3 ~~ inc_w_t3
  inc_w_t4 ~~ inc_w_t4

  gpd_w_t1 ~~ gpd_w_t1
  gpd_w_t2 ~~ gpd_w_t2
  gpd_w_t3 ~~ gpd_w_t3
  gpd_w_t4 ~~ gpd_w_t4

  # the growth factors are unrelated to the wave-one residual states
  I_inc ~~ 0*inc_w_t1
  I_inc ~~ 0*gpd_w_t1
  S_inc ~~ 0*inc_w_t1
  S_inc ~~ 0*gpd_w_t1
  I_gpd ~~ 0*inc_w_t1
  I_gpd ~~ 0*gpd_w_t1
  S_gpd ~~ 0*inc_w_t1
  S_gpd ~~ 0*gpd_w_t1

  # autoregressive effects
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

  # within-wave residual covariances
  inc_w_t1 ~~ gpd_w_t1
  inc_w_t2 ~~ gpd_w_t2
  inc_w_t3 ~~ gpd_w_t3
  inc_w_t4 ~~ gpd_w_t4
'

# Constrained: autoregressive and cross-lagged paths held equal across
# waves, using the same a/b/c/d labels as the RI-CLPM.
model_lcmsr_4w_constrained <- '

  # random-intercept and random-slope factors
  I_inc =~ 1*ink1 + 1*ink2 + 1*ink3 + 1*ink4
  S_inc =~ 0*ink1 + 1*ink2 + 2*ink3 + 3*ink4

  I_gpd =~ 1*gpd1 + 1*gpd2 + 1*gpd3 + 1*gpd4
  S_gpd =~ 0*gpd1 + 1*gpd2 + 2*gpd3 + 3*gpd4

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
  ink1 ~ 0*1
  ink2 ~ 0*1
  ink3 ~ 0*1
  ink4 ~ 0*1

  gpd1 ~~ 0*gpd1
  gpd2 ~~ 0*gpd2
  gpd3 ~~ 0*gpd3
  gpd4 ~~ 0*gpd4
  gpd1 ~ 0*1
  gpd2 ~ 0*1
  gpd3 ~ 0*1
  gpd4 ~ 0*1

  inc_w_t1 ~ 0*1
  inc_w_t2 ~ 0*1
  inc_w_t3 ~ 0*1
  inc_w_t4 ~ 0*1

  gpd_w_t1 ~ 0*1
  gpd_w_t2 ~ 0*1
  gpd_w_t3 ~ 0*1
  gpd_w_t4 ~ 0*1

  I_inc ~ 1
  S_inc ~ 1
  I_gpd ~ 1
  S_gpd ~ 1

  I_inc ~~ I_inc
  S_inc ~~ S_inc
  I_gpd ~~ I_gpd
  S_gpd ~~ S_gpd

  I_inc ~~ S_inc
  I_inc ~~ I_gpd
  I_inc ~~ S_gpd
  S_inc ~~ I_gpd
  S_inc ~~ S_gpd
  I_gpd ~~ S_gpd

  inc_w_t1 ~~ inc_w_t1
  inc_w_t2 ~~ inc_w_t2
  inc_w_t3 ~~ inc_w_t3
  inc_w_t4 ~~ inc_w_t4

  gpd_w_t1 ~~ gpd_w_t1
  gpd_w_t2 ~~ gpd_w_t2
  gpd_w_t3 ~~ gpd_w_t3
  gpd_w_t4 ~~ gpd_w_t4

  I_inc ~~ 0*inc_w_t1
  I_inc ~~ 0*gpd_w_t1
  S_inc ~~ 0*inc_w_t1
  S_inc ~~ 0*gpd_w_t1
  I_gpd ~~ 0*inc_w_t1
  I_gpd ~~ 0*gpd_w_t1
  S_gpd ~~ 0*inc_w_t1
  S_gpd ~~ 0*gpd_w_t1

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
'

########## Define the LCM-SR: 6 waves (Trier) ####################################

model_lcmsr_6w <- '

  # random-intercept and random-slope factors
  I_inc =~ 1*ink1 + 1*ink2 + 1*ink3 + 1*ink4 + 1*ink5 + 1*ink6
  S_inc =~ 0*ink1 + 1*ink2 + 2*ink3 + 3*ink4 + 4*ink5 + 5*ink6

  I_gpd =~ 1*gpd1 + 1*gpd2 + 1*gpd3 + 1*gpd4 + 1*gpd5 + 1*gpd6
  S_gpd =~ 0*gpd1 + 1*gpd2 + 2*gpd3 + 3*gpd4 + 4*gpd5 + 5*gpd6

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
  ink1 ~ 0*1
  ink2 ~ 0*1
  ink3 ~ 0*1
  ink4 ~ 0*1
  ink5 ~ 0*1
  ink6 ~ 0*1

  gpd1 ~~ 0*gpd1
  gpd2 ~~ 0*gpd2
  gpd3 ~~ 0*gpd3
  gpd4 ~~ 0*gpd4
  gpd5 ~~ 0*gpd5
  gpd6 ~~ 0*gpd6
  gpd1 ~ 0*1
  gpd2 ~ 0*1
  gpd3 ~ 0*1
  gpd4 ~ 0*1
  gpd5 ~ 0*1
  gpd6 ~ 0*1

  inc_w_t1 ~ 0*1
  inc_w_t2 ~ 0*1
  inc_w_t3 ~ 0*1
  inc_w_t4 ~ 0*1
  inc_w_t5 ~ 0*1
  inc_w_t6 ~ 0*1

  gpd_w_t1 ~ 0*1
  gpd_w_t2 ~ 0*1
  gpd_w_t3 ~ 0*1
  gpd_w_t4 ~ 0*1
  gpd_w_t5 ~ 0*1
  gpd_w_t6 ~ 0*1

  I_inc ~ 1
  S_inc ~ 1
  I_gpd ~ 1
  S_gpd ~ 1

  I_inc ~~ I_inc
  S_inc ~~ S_inc
  I_gpd ~~ I_gpd
  S_gpd ~~ S_gpd

  I_inc ~~ S_inc
  I_inc ~~ I_gpd
  I_inc ~~ S_gpd
  S_inc ~~ I_gpd
  S_inc ~~ S_gpd
  I_gpd ~~ S_gpd

  inc_w_t1 ~~ inc_w_t1
  inc_w_t2 ~~ inc_w_t2
  inc_w_t3 ~~ inc_w_t3
  inc_w_t4 ~~ inc_w_t4
  inc_w_t5 ~~ inc_w_t5
  inc_w_t6 ~~ inc_w_t6

  gpd_w_t1 ~~ gpd_w_t1
  gpd_w_t2 ~~ gpd_w_t2
  gpd_w_t3 ~~ gpd_w_t3
  gpd_w_t4 ~~ gpd_w_t4
  gpd_w_t5 ~~ gpd_w_t5
  gpd_w_t6 ~~ gpd_w_t6

  I_inc ~~ 0*inc_w_t1
  I_inc ~~ 0*gpd_w_t1
  S_inc ~~ 0*inc_w_t1
  S_inc ~~ 0*gpd_w_t1
  I_gpd ~~ 0*inc_w_t1
  I_gpd ~~ 0*gpd_w_t1
  S_gpd ~~ 0*inc_w_t1
  S_gpd ~~ 0*gpd_w_t1

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
'

model_lcmsr_6w_constrained <- '

  # random-intercept and random-slope factors
  I_inc =~ 1*ink1 + 1*ink2 + 1*ink3 + 1*ink4 + 1*ink5 + 1*ink6
  S_inc =~ 0*ink1 + 1*ink2 + 2*ink3 + 3*ink4 + 4*ink5 + 5*ink6

  I_gpd =~ 1*gpd1 + 1*gpd2 + 1*gpd3 + 1*gpd4 + 1*gpd5 + 1*gpd6
  S_gpd =~ 0*gpd1 + 1*gpd2 + 2*gpd3 + 3*gpd4 + 4*gpd5 + 5*gpd6

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
  ink1 ~ 0*1
  ink2 ~ 0*1
  ink3 ~ 0*1
  ink4 ~ 0*1
  ink5 ~ 0*1
  ink6 ~ 0*1

  gpd1 ~~ 0*gpd1
  gpd2 ~~ 0*gpd2
  gpd3 ~~ 0*gpd3
  gpd4 ~~ 0*gpd4
  gpd5 ~~ 0*gpd5
  gpd6 ~~ 0*gpd6
  gpd1 ~ 0*1
  gpd2 ~ 0*1
  gpd3 ~ 0*1
  gpd4 ~ 0*1
  gpd5 ~ 0*1
  gpd6 ~ 0*1

  inc_w_t1 ~ 0*1
  inc_w_t2 ~ 0*1
  inc_w_t3 ~ 0*1
  inc_w_t4 ~ 0*1
  inc_w_t5 ~ 0*1
  inc_w_t6 ~ 0*1

  gpd_w_t1 ~ 0*1
  gpd_w_t2 ~ 0*1
  gpd_w_t3 ~ 0*1
  gpd_w_t4 ~ 0*1
  gpd_w_t5 ~ 0*1
  gpd_w_t6 ~ 0*1

  I_inc ~ 1
  S_inc ~ 1
  I_gpd ~ 1
  S_gpd ~ 1

  I_inc ~~ I_inc
  S_inc ~~ S_inc
  I_gpd ~~ I_gpd
  S_gpd ~~ S_gpd

  I_inc ~~ S_inc
  I_inc ~~ I_gpd
  I_inc ~~ S_gpd
  S_inc ~~ I_gpd
  S_inc ~~ S_gpd
  I_gpd ~~ S_gpd

  inc_w_t1 ~~ inc_w_t1
  inc_w_t2 ~~ inc_w_t2
  inc_w_t3 ~~ inc_w_t3
  inc_w_t4 ~~ inc_w_t4
  inc_w_t5 ~~ inc_w_t5
  inc_w_t6 ~~ inc_w_t6

  gpd_w_t1 ~~ gpd_w_t1
  gpd_w_t2 ~~ gpd_w_t2
  gpd_w_t3 ~~ gpd_w_t3
  gpd_w_t4 ~~ gpd_w_t4
  gpd_w_t5 ~~ gpd_w_t5
  gpd_w_t6 ~~ gpd_w_t6

  I_inc ~~ 0*inc_w_t1
  I_inc ~~ 0*gpd_w_t1
  S_inc ~~ 0*inc_w_t1
  S_inc ~~ 0*gpd_w_t1
  I_gpd ~~ 0*inc_w_t1
  I_gpd ~~ 0*gpd_w_t1
  S_gpd ~~ 0*inc_w_t1
  S_gpd ~~ 0*gpd_w_t1

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
'

########## Standardized effects with 95% CIs ###################################

# Same helpers as in 03_riclpm.R, with the LCM-SR's intercept factors (I_inc,
# I_gpd) in place of the RI-CLPM's random intercepts.
std_ci_table <- function(fit, intercepts = c("I_inc", "I_gpd")) {
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

# The values reported in Table S4: each autoregressive and cross-lagged effect,
# and its confidence limits, averaged across waves (the constrained paths
# differ across waves only in their standardization), plus the correlation
# between the intercept factors.
std_ci_averaged <- function(fit, intercepts = c("I_inc", "I_gpd")) {
  paths <- std_ci_table(fit, intercepts)

  paths$effect <- factor(
    ifelse(paths$op == "~~", "r_INC,GPD",
    ifelse(grepl("^gpd_w", paths$lhs) & grepl("^inc_w", paths$rhs), "INC -> GPD",
    ifelse(grepl("^inc_w", paths$lhs) & grepl("^gpd_w", paths$rhs), "GPD -> INC",
    ifelse(grepl("^inc_w", paths$lhs), "INC -> INC", "GPD -> GPD")))),
    levels = c("r_INC,GPD", "INC -> GPD", "GPD -> INC", "INC -> INC", "GPD -> GPD"))

  out <- paths %>%
    group_by(effect) %>%
    summarise(est.std = mean(est.std),
              ci.lower = mean(ci.lower),
              ci.upper = mean(ci.upper),
              .groups = "drop") %>%
    as.data.frame()
  out[, -1] <- round(out[, -1], 6)
  out
}

print_std_ci <- function(fit, intercepts = c("I_inc", "I_gpd"), averaged = FALSE) {
  cat("\n\nStandardized (std.all) regression paths and intercept-factor\n")
  cat("correlation (", intercepts[1], " ~~ ", intercepts[2], "), with 95% ",
      "confidence intervals:\n\n", sep = "")
  if (!lavInspect(fit, "converged")) {
    cat("WARNING: this model did not converge, so the values below are ",
        "unreliable (NA = not computable).\n\n", sep = "")
  }
  print(std_ci_table(fit, intercepts))

  if (averaged) {
    cat("\n\nThe same effects averaged across waves (estimates and confidence\n")
    cat("limits), as reported in Table S4:\n\n")
    print(std_ci_averaged(fit, intercepts), row.names = FALSE)
  }
}

########## Bern: unconstrained model ##############################################

fit_lcmsr_bern <- sem(
  model_lcmsr_4w,
  data = df_bern,
  missing = "ML")

summary(fit_lcmsr_bern, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_lcmsr_bern)

lavInspect(fit_lcmsr_bern, "post.check")

########## Bern: constrained model #################################################

fit_lcmsr_bern_constrained <- sem(
  model_lcmsr_4w_constrained,
  data = df_bern,
  missing = "ML")

summary(fit_lcmsr_bern_constrained, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_lcmsr_bern_constrained)

lavInspect(fit_lcmsr_bern_constrained, "post.check")

########## Bern: constrained vs. unconstrained #####################################

if (lavInspect(fit_lcmsr_bern, "converged")) {
  summary(compareFit(fit_lcmsr_bern_constrained, fit_lcmsr_bern))

  fitMeasures(
    fit_lcmsr_bern_constrained,
    "rmsea",
    h1.model = fit_lcmsr_bern)
} else {
  message("Bern: the unconstrained LCM-SR did not converge, so the comparison is skipped.")
}

########## Zurich: unconstrained model #############################################

fit_lcmsr_zurich <- sem(
  model_lcmsr_4w,
  data = df_zurich,
  missing = "ML")

summary(fit_lcmsr_zurich, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_lcmsr_zurich)

lavInspect(fit_lcmsr_zurich, "post.check")

########## Zurich: constrained model ################################################

fit_lcmsr_zurich_constrained <- sem(
  model_lcmsr_4w_constrained,
  data = df_zurich,
  missing = "ML")

summary(fit_lcmsr_zurich_constrained, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_lcmsr_zurich_constrained)

lavInspect(fit_lcmsr_zurich_constrained, "post.check")

########## Zurich: constrained vs. unconstrained ####################################

if (lavInspect(fit_lcmsr_zurich, "converged")) {
  summary(compareFit(fit_lcmsr_zurich_constrained, fit_lcmsr_zurich))

  fitMeasures(
    fit_lcmsr_zurich_constrained,
    "rmsea",
    h1.model = fit_lcmsr_zurich)
} else {
  message("Zurich: the unconstrained LCM-SR did not converge, so the comparison is skipped.")
}

########## Trier: unconstrained model, 6 waves #######################################

fit_lcmsr_trier <- sem(
  model_lcmsr_6w,
  data = df_trier,
  missing = "ML")

summary(fit_lcmsr_trier, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_lcmsr_trier)

lavInspect(fit_lcmsr_trier, "post.check")

########## Trier: constrained model, 6 waves #########################################

fit_lcmsr_trier_constrained <- sem(
  model_lcmsr_6w_constrained,
  data = df_trier,
  missing = "ML")

summary(fit_lcmsr_trier_constrained, standardized = TRUE, fit.measures = TRUE)

std_ci_table(fit_lcmsr_trier_constrained)

lavInspect(fit_lcmsr_trier_constrained, "post.check")

########## Trier: constrained vs. unconstrained #######################################

if (lavInspect(fit_lcmsr_trier, "converged")) {
  summary(compareFit(fit_lcmsr_trier_constrained, fit_lcmsr_trier))

  fitMeasures(
    fit_lcmsr_trier_constrained,
    "rmsea",
    h1.model = fit_lcmsr_trier)
} else {
  message("Trier: the unconstrained LCM-SR did not converge, so the comparison is skipped.")
}

########## Admissibility check ########################################################

# Whether each model reached a proper (admissible) solution: lavaan's own
# post-estimation check, and whether any estimated latent variance came out
# negative.

lcmsr_admissibility_check <- tibble(
  sample = c("Bern", "Bern", "Zurich", "Zurich", "Trier", "Trier"),
  model = c(
    "constrained", "unconstrained",
    "constrained", "unconstrained",
    "constrained", "unconstrained"),
  post_check = c(
    lavInspect(fit_lcmsr_bern_constrained, "post.check"),
    lavInspect(fit_lcmsr_bern, "post.check"),
    lavInspect(fit_lcmsr_zurich_constrained, "post.check"),
    lavInspect(fit_lcmsr_zurich, "post.check"),
    lavInspect(fit_lcmsr_trier_constrained, "post.check"),
    lavInspect(fit_lcmsr_trier, "post.check")),
  negative_variances = c(
    sum(diag(lavInspect(fit_lcmsr_bern_constrained, "cov.lv")) < 0),
    sum(diag(lavInspect(fit_lcmsr_bern, "cov.lv")) < 0),
    sum(diag(lavInspect(fit_lcmsr_zurich_constrained, "cov.lv")) < 0),
    sum(diag(lavInspect(fit_lcmsr_zurich, "cov.lv")) < 0),
    sum(diag(lavInspect(fit_lcmsr_trier_constrained, "cov.lv")) < 0),
    sum(diag(lavInspect(fit_lcmsr_trier, "cov.lv")) < 0))
)

print(lcmsr_admissibility_check)

if (any(lcmsr_admissibility_check$post_check == FALSE)) {
  message(
    "Note: at least one LCM-SR is not admissible (post.check == FALSE / ",
    "a negative latent variance). It is still reported above, as a ",
    "transparent exploratory result -- this is the reason the LCM-SR is ",
    "kept in the supplement rather than treated as the main model.")
}

################################################################################
########## Save model output as text files ####################################
#
# Each model's raw R/lavaan output (the unconstrained summary, the
# constrained summary, and the constrained-vs-unconstrained comparison),
# exactly as summary()/compareFit()/fitMeasures() produce it above, plus the
# admissibility check. Each summary file ends with the standardized paths' 95%
# CI table; the constrained summaries also end with the wave-averaged values
# reported in Table S4.

summaries_dir <- "results/04 lcm-sr exploratory analysis/model summaries"
comparisons_dir <- "results/04 lcm-sr exploratory analysis/model comparisons"
dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparisons_dir, recursive = TRUE, showWarnings = FALSE)

save_output <- function(file, expr) {
  writeLines(capture.output(force(expr)), con = file)
}

# The unconstrained LCM-SR does not always converge; compareFit() and
# fitMeasures(h1.model = ) cannot be trusted on a model that did not
# converge, so the comparison is written only when both models converged,
# and a plain explanation is written otherwise.
save_comparison <- function(file, fit_constrained, fit_unconstrained, sample_name) {
  if (lavInspect(fit_unconstrained, "converged")) {
    save_output(file, {
      print(summary(compareFit(fit_constrained, fit_unconstrained)))
      cat("\n\nRMSEA using unconstrained model as h1 model:\n")
      print(fitMeasures(fit_constrained, "rmsea", h1.model = fit_unconstrained))
    })
  } else {
    save_output(file, {
      cat(
        sample_name, ": the unconstrained model did not converge, so the ",
        "constrained vs. unconstrained comparison could not be computed. ",
        "See the constrained model's own summary file instead.\n",
        sep = ""
      )
    })
  }
}

save_output(
  file.path(summaries_dir, "lcm_sr_bern_4w_unconstrained_summary.txt"),
  {
    print(summary(fit_lcmsr_bern, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_lcmsr_bern)
  }
)
save_output(
  file.path(summaries_dir, "lcm_sr_bern_4w_constrained_summary.txt"),
  {
    print(summary(fit_lcmsr_bern_constrained, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_lcmsr_bern_constrained, averaged = TRUE)
  }
)
save_comparison(
  file.path(comparisons_dir, "lcm_sr_bern_4w_model_comparison.txt"),
  fit_lcmsr_bern_constrained, fit_lcmsr_bern, "Bern"
)

save_output(
  file.path(summaries_dir, "lcm_sr_zurich_4w_unconstrained_summary.txt"),
  {
    print(summary(fit_lcmsr_zurich, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_lcmsr_zurich)
  }
)
save_output(
  file.path(summaries_dir, "lcm_sr_zurich_4w_constrained_summary.txt"),
  {
    print(summary(fit_lcmsr_zurich_constrained, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_lcmsr_zurich_constrained, averaged = TRUE)
  }
)
save_comparison(
  file.path(comparisons_dir, "lcm_sr_zurich_4w_model_comparison.txt"),
  fit_lcmsr_zurich_constrained, fit_lcmsr_zurich, "Zurich"
)

save_output(
  file.path(summaries_dir, "lcm_sr_trier_6w_unconstrained_summary.txt"),
  {
    print(summary(fit_lcmsr_trier, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_lcmsr_trier)
  }
)
save_output(
  file.path(summaries_dir, "lcm_sr_trier_6w_constrained_summary.txt"),
  {
    print(summary(fit_lcmsr_trier_constrained, standardized = TRUE, fit.measures = TRUE))
    print_std_ci(fit_lcmsr_trier_constrained, averaged = TRUE)
  }
)
save_comparison(
  file.path(comparisons_dir, "lcm_sr_trier_6w_model_comparison.txt"),
  fit_lcmsr_trier_constrained, fit_lcmsr_trier, "Trier"
)

save_output(
  file.path(comparisons_dir, "lcm_sr_admissibility_check.txt"),
  print(lcmsr_admissibility_check)
)

################################################################################
########## Plot the constrained LCM-SRs #######################################

# Same idea as the RI-CLPM layout in the main analysis, but with two growth
# factors per construct (S above, I below) instead of one.
lay_lcmsr_4w <- get_layout(
  "S_inc",    NA,         NA,         "I_inc",
  "ink1",     "ink2",     "ink3",     "ink4",
  "inc_w_t1", "inc_w_t2", "inc_w_t3", "inc_w_t4",
  "gpd_w_t1", "gpd_w_t2", "gpd_w_t3", "gpd_w_t4",
  "gpd1",     "gpd2",     "gpd3",     "gpd4",
  "S_gpd",    NA,         NA,         "I_gpd",
  rows = 6)

lay_lcmsr_6w <- get_layout(
  "S_inc",    NA,         NA,         NA,         NA,         "I_inc",
  "ink1",     "ink2",     "ink3",     "ink4",     "ink5",     "ink6",
  "inc_w_t1", "inc_w_t2", "inc_w_t3", "inc_w_t4", "inc_w_t5", "inc_w_t6",
  "gpd_w_t1", "gpd_w_t2", "gpd_w_t3", "gpd_w_t4", "gpd_w_t5", "gpd_w_t6",
  "gpd1",     "gpd2",     "gpd3",     "gpd4",     "gpd5",     "gpd6",
  "S_gpd",    NA,         NA,         NA,         NA,         "I_gpd",
  rows = 6)

########## Bern plots #################################################################

p_lcmsr_bern_names <- graph_sem(
  fit_lcmsr_bern_constrained,
  layout = lay_lcmsr_4w,
  label = "name")

ggsave(
  "plots/04 lcm-sr exploratory analysis/lcm_sr_bern_names.png",
  plot = p_lcmsr_bern_names,
  width = 13, height = 10, dpi = 300, bg = "white")

p_lcmsr_bern_est <- graph_sem(
  fit_lcmsr_bern_constrained,
  layout = lay_lcmsr_4w,
  label = "est_std")

ggsave(
  "plots/04 lcm-sr exploratory analysis/lcm_sr_bern_est_std.png",
  plot = p_lcmsr_bern_est,
  width = 13, height = 10, dpi = 300, bg = "white")

########## Zurich plots #################################################################

p_lcmsr_zurich_names <- graph_sem(
  fit_lcmsr_zurich_constrained,
  layout = lay_lcmsr_4w,
  label = "name")

ggsave(
  "plots/04 lcm-sr exploratory analysis/lcm_sr_zurich_names.png",
  plot = p_lcmsr_zurich_names,
  width = 13, height = 10, dpi = 300, bg = "white")

p_lcmsr_zurich_est <- graph_sem(
  fit_lcmsr_zurich_constrained,
  layout = lay_lcmsr_4w,
  label = "est_std")

ggsave(
  "plots/04 lcm-sr exploratory analysis/lcm_sr_zurich_est_std.png",
  plot = p_lcmsr_zurich_est,
  width = 13, height = 10, dpi = 300, bg = "white")

########## Trier plots ####################################################################

p_lcmsr_trier_names <- graph_sem(
  fit_lcmsr_trier_constrained,
  layout = lay_lcmsr_6w,
  label = "name")

ggsave(
  "plots/04 lcm-sr exploratory analysis/lcm_sr_trier_names.png",
  plot = p_lcmsr_trier_names,
  width = 15, height = 10, dpi = 300, bg = "white")

p_lcmsr_trier_est <- graph_sem(
  fit_lcmsr_trier_constrained,
  layout = lay_lcmsr_6w,
  label = "est_std")

ggsave(
  "plots/04 lcm-sr exploratory analysis/lcm_sr_trier_est_std.png",
  plot = p_lcmsr_trier_est,
  width = 15, height = 10, dpi = 300, bg = "white")
