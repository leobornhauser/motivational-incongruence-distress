# Calculates means, SDs, and available sample sizes for INC and GPD.
# Run from the project root directory after building the analysis datasets.

library(readr)
library(dplyr)
library(tidyr)

########## Data ####################################################################

bern_scores <- read_csv(
  "data/processed/bern_scores.csv",
  col_types = cols(patient_id = col_character()),
  show_col_types = FALSE
)

zurich_scores <- read_csv(
  "data/processed/zurich_scores.csv",
  col_types = cols(patient_id = col_character()),
  show_col_types = FALSE
)

trier_scores <- read_csv(
  "data/processed/trier_scores.csv",
  col_types = cols(patient_id = col_character()),
  show_col_types = FALSE
)

########## Patient identifiers #####################################################

score_ids <- bind_rows(
  bern_scores %>% transmute(sample = "Bern", patient_id),
  zurich_scores %>% transmute(sample = "Zurich", patient_id),
  trier_scores %>% transmute(sample = "Trier", patient_id)
)

stopifnot(
  !anyNA(score_ids$patient_id),
  all(nzchar(trimws(score_ids$patient_id))),
  anyDuplicated(score_ids[c("sample", "patient_id")]) == 0L
)

########## Score variables #########################################################

trial_variables <- c(
  "ink_total_pre", "ink_total_z1", "ink_total_z2", "ink_total_post",
  "gpd_total_pre", "gpd_total_z1", "gpd_total_z2", "gpd_total_post"
)

trier_variables <- c(
  "ink_total_pre", "ink_total_s05", "ink_total_s10",
  "ink_total_s15", "ink_total_s20", "ink_total_s25",
  "gpd_total_pre", "gpd_total_s05", "gpd_total_s10",
  "gpd_total_s15", "gpd_total_s20", "gpd_total_s25"
)

########## Long-format scores ######################################################

bern_long <- bern_scores %>%
  select(all_of(trial_variables)) %>%
  pivot_longer(
    everything(),
    names_to = c("construct", "wave"),
    names_pattern = "(ink|gpd)_total_(.*)",
    values_to = "score"
  ) %>%
  mutate(sample = "Bern")

zurich_long <- zurich_scores %>%
  select(all_of(trial_variables)) %>%
  pivot_longer(
    everything(),
    names_to = c("construct", "wave"),
    names_pattern = "(ink|gpd)_total_(.*)",
    values_to = "score"
  ) %>%
  mutate(sample = "Zurich")

trier_long <- trier_scores %>%
  select(all_of(trier_variables)) %>%
  pivot_longer(
    everything(),
    names_to = c("construct", "wave"),
    names_pattern = "(ink|gpd)_total_(.*)",
    values_to = "score"
  ) %>%
  mutate(sample = "Trier")

wave_scores <- bind_rows(
  bern_long,
  zurich_long,
  trier_long
) %>%
  mutate(
    construct = recode(construct, ink = "INC", gpd = "GPD"),
    
    wave_order = case_when(
      wave == "pre" ~ 1L,
      wave %in% c("z1", "s05") ~ 2L,
      wave %in% c("z2", "s10") ~ 3L,
      wave %in% c("post", "s15") ~ 4L,
      wave == "s20" ~ 5L,
      wave == "s25" ~ 6L
    ),
    
    assessment = case_when(
      wave == "pre" ~ "Pretreatment",
      sample == "Bern" & wave == "z1" ~ "Session 8",
      sample == "Bern" & wave == "z2" ~ "Session 16",
      sample == "Zurich" & wave == "z1" ~ "Session 7",
      sample == "Zurich" & wave == "z2" ~ "Session 14",
      wave == "post" ~ "Posttreatment",
      wave == "s05" ~ "Session 5",
      wave == "s10" ~ "Session 10",
      wave == "s15" ~ "Session 15",
      wave == "s20" ~ "Session 20",
      wave == "s25" ~ "Session 25"
    ),
    
    sample = factor(sample, levels = c("Bern", "Zurich", "Trier"))
  )

########## Wave descriptives #######################################################

# Each mean and SD uses the available scores for that construct and wave.
wave_descriptives <- wave_scores %>%
  group_by(sample, wave_order, assessment, construct) %>%
  summarise(
    analytic_n = n(),
    available_n = sum(!is.na(score)),
    missing_n = sum(is.na(score)),
    
    mean = if (any(!is.na(score))) {
      mean(score, na.rm = TRUE)
    } else {
      NA_real_
    },
    
    sd = if (sum(!is.na(score)) > 1L) {
      sd(score, na.rm = TRUE)
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  arrange(sample, wave_order, construct)

########## Table 1 #################################################################

# Format display values only; retain full precision in wave_descriptives.
table_1 <- wave_descriptives %>%
  mutate(
    mean_sd = case_when(
      available_n == 0L ~ "NA",
      available_n == 1L ~ sprintf("%.2f (NA)", mean),
      TRUE ~ sprintf("%.2f (%.2f)", mean, sd)
    )
  ) %>%
  select(
    sample,
    analytic_n,
    wave_order,
    assessment,
    construct,
    mean_sd,
    available_n
  ) %>%
  pivot_wider(
    names_from = construct,
    values_from = c(mean_sd, available_n),
    names_glue = "{construct}_{.value}"
  ) %>%
  arrange(sample, wave_order) %>%
  select(
    sample,
    analytic_n,
    assessment,
    INC_mean_sd,
    INC_n = INC_available_n,
    GPD_mean_sd,
    GPD_n = GPD_available_n
  )

########## Output ##################################################################

print(table_1, n = Inf, width = Inf)

print(
  wave_descriptives %>%
    select(
      sample,
      assessment,
      construct,
      analytic_n,
      available_n,
      missing_n
    ),
  n = Inf,
  width = Inf
)