# Calculates raw Cronbach's alpha using complete baseline item sets.
# Items are already prepared; no further reversal or transformation is applied.
# Run from the project root directory after building the analysis datasets.

library(readr)
library(dplyr)
library(psych)

########## Prepared datasets #######################################################

bern_items <- read_csv(
  "data_raw/bern_master_inc_gpd.csv",
  col_types = cols(patient_id = col_character()),
  show_col_types = FALSE
)

zurich_items <- read_csv(
  "data_raw/zurich_master_inc_gpd.csv",
  col_types = cols(patient_id = col_character()),
  show_col_types = FALSE
)

trier_items <- read_csv(
  "data_raw/trier_master_inc_gpd.csv",
  col_types = cols(patient_id = col_character()),
  show_col_types = FALSE
)

########## Analytic score datasets #################################################

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

input_ids <- bind_rows(
  bern_items %>% transmute(source = "Bern items", patient_id),
  bern_scores %>% transmute(source = "Bern scores", patient_id),
  zurich_items %>% transmute(source = "Zurich items", patient_id),
  zurich_scores %>% transmute(source = "Zurich scores", patient_id),
  trier_items %>% transmute(source = "Trier items", patient_id),
  trier_scores %>% transmute(source = "Trier scores", patient_id)
)

stopifnot(
  !anyNA(input_ids$patient_id),
  all(nzchar(trimws(input_ids$patient_id))),
  anyDuplicated(input_ids[c("source", "patient_id")]) == 0L,
  all(bern_scores$patient_id %in% bern_items$patient_id),
  all(zurich_scores$patient_id %in% zurich_items$patient_id),
  all(trier_scores$patient_id %in% trier_items$patient_id)
)

########## Analytic samples ########################################################

bern_analysis <- bern_scores %>%
  select(patient_id) %>%
  left_join(bern_items, by = "patient_id")

zurich_analysis <- zurich_scores %>%
  select(patient_id) %>%
  left_join(zurich_items, by = "patient_id")

trier_analysis <- trier_scores %>%
  select(patient_id) %>%
  left_join(trier_items, by = "patient_id")

########## Baseline items ##########################################################

bern_inc <- bern_analysis %>%
  select(all_of(sprintf("inc_%02d_pre", 1:23)))

bern_gpd <- bern_analysis %>%
  select(all_of(sprintf("bsi_%02d_pre", 1:53)))

zurich_inc <- zurich_analysis %>%
  select(all_of(sprintf("inc_%02d_pre", 1:23)))

zurich_gpd <- zurich_analysis %>%
  select(all_of(sprintf("scl_%02d_pre", 1:9)))

trier_inc <- trier_analysis %>%
  select(all_of(sprintf("inc_%02d_pre", 1:10)))

# Baseline GPD uses eleven mapped BSI items on the 1-4 scale.
trier_gpd <- trier_analysis %>%
  select(all_of(sprintf("gpd_%02d_pre", 1:11)))

########## Complete item sets ######################################################

# Select complete cases separately for each instrument.
bern_inc_complete <- bern_inc[complete.cases(bern_inc), ]
bern_gpd_complete <- bern_gpd[complete.cases(bern_gpd), ]

zurich_inc_complete <- zurich_inc[complete.cases(zurich_inc), ]
zurich_gpd_complete <- zurich_gpd[complete.cases(zurich_gpd), ]

trier_inc_complete <- trier_inc[complete.cases(trier_inc), ]
trier_gpd_complete <- trier_gpd[complete.cases(trier_gpd), ]

# Item completeness must agree with baseline score availability.
stopifnot(
  all(complete.cases(bern_inc) == !is.na(bern_scores$ink_total_pre)),
  all(complete.cases(bern_gpd) == !is.na(bern_scores$gpd_total_pre)),
  all(complete.cases(zurich_inc) == !is.na(zurich_scores$ink_total_pre)),
  all(complete.cases(zurich_gpd) == !is.na(zurich_scores$gpd_total_pre)),
  all(complete.cases(trier_inc) == !is.na(trier_scores$ink_total_pre)),
  all(complete.cases(trier_gpd) == !is.na(trier_scores$gpd_total_pre))
)

# Require at least two observations and nonconstant items.
stopifnot(
  nrow(bern_inc_complete) > 1L,
  nrow(bern_gpd_complete) > 1L,
  nrow(zurich_inc_complete) > 1L,
  nrow(zurich_gpd_complete) > 1L,
  nrow(trier_inc_complete) > 1L,
  nrow(trier_gpd_complete) > 1L,
  all(vapply(bern_inc_complete, var, numeric(1)) > 0),
  all(vapply(bern_gpd_complete, var, numeric(1)) > 0),
  all(vapply(zurich_inc_complete, var, numeric(1)) > 0),
  all(vapply(zurich_gpd_complete, var, numeric(1)) > 0),
  all(vapply(trier_inc_complete, var, numeric(1)) > 0),
  all(vapply(trier_gpd_complete, var, numeric(1)) > 0)
)

########## Bern trial ##############################################################

# Calculate ordinary item alpha without automatic reversal or item removal.
alpha_bern_inc <- psych::alpha(
  bern_inc_complete,
  check.keys = FALSE,
  delete = FALSE,
  use = "all.obs"
)

alpha_bern_gpd <- psych::alpha(
  bern_gpd_complete,
  check.keys = FALSE,
  delete = FALSE,
  use = "all.obs"
)

########## Zurich trial ############################################################

alpha_zurich_inc <- psych::alpha(
  zurich_inc_complete,
  check.keys = FALSE,
  delete = FALSE,
  use = "all.obs"
)

alpha_zurich_gpd <- psych::alpha(
  zurich_gpd_complete,
  check.keys = FALSE,
  delete = FALSE,
  use = "all.obs"
)

########## Trier sample ############################################################

alpha_trier_inc <- psych::alpha(
  trier_inc_complete,
  check.keys = FALSE,
  delete = FALSE,
  use = "all.obs"
)

alpha_trier_gpd <- psych::alpha(
  trier_gpd_complete,
  check.keys = FALSE,
  delete = FALSE,
  use = "all.obs",
  discrete = FALSE
)

########## Reliability summary #####################################################

baseline_reliability <- tibble(
  sample = c("Bern", "Bern", "Zurich", "Zurich", "Trier", "Trier"),
  construct = c("INC", "GPD", "INC", "GPD", "INC", "GPD"),
  
  item_set = c(
    "INK-23",
    "BSI-53",
    "INK-23",
    "SCL-K-9",
    "INK-10",
    "11 mapped baseline BSI items"
  ),
  
  analytic_n = c(
    nrow(bern_analysis), nrow(bern_analysis),
    nrow(zurich_analysis), nrow(zurich_analysis),
    nrow(trier_analysis), nrow(trier_analysis)
  ),
  
  n_items = c(
    ncol(bern_inc), ncol(bern_gpd),
    ncol(zurich_inc), ncol(zurich_gpd),
    ncol(trier_inc), ncol(trier_gpd)
  ),
  
  alpha_n = c(
    nrow(bern_inc_complete), nrow(bern_gpd_complete),
    nrow(zurich_inc_complete), nrow(zurich_gpd_complete),
    nrow(trier_inc_complete), nrow(trier_gpd_complete)
  ),
  
  raw_alpha = c(
    alpha_bern_inc$total$raw_alpha,
    alpha_bern_gpd$total$raw_alpha,
    alpha_zurich_inc$total$raw_alpha,
    alpha_zurich_gpd$total$raw_alpha,
    alpha_trier_inc$total$raw_alpha,
    alpha_trier_gpd$total$raw_alpha
  )
) %>%
  mutate(n_without_complete_items = analytic_n - alpha_n)

stopifnot(all(is.finite(baseline_reliability$raw_alpha)))

########## Output ##################################################################

print(
  baseline_reliability %>%
    mutate(raw_alpha = sprintf("%.2f", raw_alpha)),
  n = Inf,
  width = Inf
)