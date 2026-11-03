# Builds analytic total-score datasets from the prepared item datasets.
# Run from the project root directory.

library(readr)
library(dplyr)

########## Data ####################################################################

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

########## Patient identifiers #####################################################

input_ids <- bind_rows(
  bern_items %>% transmute(sample = "Bern", patient_id),
  zurich_items %>% transmute(sample = "Zurich", patient_id),
  trier_items %>% transmute(sample = "Trier", patient_id)
)

stopifnot(
  !anyNA(input_ids$patient_id),
  all(nzchar(trimws(input_ids$patient_id))),
  anyDuplicated(input_ids[c("sample", "patient_id")]) == 0L
)

########## Scoring #################################################################

# Items have already been reversed or transformed where required.
# Each score requires complete responses to all contributing items.
# Retain patients with at least one observed score across modeled assessments.

########## Bern trial ##############################################################

# INC gives equal weight to the approach and avoidance subscale means.
# GPD is the mean of the 53 BSI items.
bern_scores <- bern_items %>%
  transmute(
    patient_id,
    
    ink_total_pre = (
      rowMeans(across(all_of(sprintf("inc_%02d_pre", 1:14))), na.rm = FALSE) +
        rowMeans(across(all_of(sprintf("inc_%02d_pre", 15:23))), na.rm = FALSE)
    ) / 2,
    
    ink_total_z1 = (
      rowMeans(across(all_of(sprintf("inc_%02d_z1", 1:14))), na.rm = FALSE) +
        rowMeans(across(all_of(sprintf("inc_%02d_z1", 15:23))), na.rm = FALSE)
    ) / 2,
    
    ink_total_z2 = (
      rowMeans(across(all_of(sprintf("inc_%02d_z2", 1:14))), na.rm = FALSE) +
        rowMeans(across(all_of(sprintf("inc_%02d_z2", 15:23))), na.rm = FALSE)
    ) / 2,
    
    ink_total_post = (
      rowMeans(across(all_of(sprintf("inc_%02d_post", 1:14))), na.rm = FALSE) +
        rowMeans(across(all_of(sprintf("inc_%02d_post", 15:23))), na.rm = FALSE)
    ) / 2,
    
    gpd_total_pre = rowMeans(
      across(all_of(sprintf("bsi_%02d_pre", 1:53))), na.rm = FALSE
    ),
    gpd_total_z1 = rowMeans(
      across(all_of(sprintf("bsi_%02d_z1", 1:53))), na.rm = FALSE
    ),
    gpd_total_z2 = rowMeans(
      across(all_of(sprintf("bsi_%02d_z2", 1:53))), na.rm = FALSE
    ),
    gpd_total_post = rowMeans(
      across(all_of(sprintf("bsi_%02d_post", 1:53))), na.rm = FALSE
    )
  ) %>%
  filter(if_any(-patient_id, ~ !is.na(.x)))

########## Zurich trial ############################################################

# INC gives equal weight to the approach and avoidance subscale means.
# GPD is the mean of the nine SCL-K-9 items.
zurich_scores <- zurich_items %>%
  transmute(
    patient_id,
    
    ink_total_pre = (
      rowMeans(across(all_of(sprintf("inc_%02d_pre", 1:14))), na.rm = FALSE) +
        rowMeans(across(all_of(sprintf("inc_%02d_pre", 15:23))), na.rm = FALSE)
    ) / 2,
    
    ink_total_z1 = (
      rowMeans(across(all_of(sprintf("inc_%02d_z1", 1:14))), na.rm = FALSE) +
        rowMeans(across(all_of(sprintf("inc_%02d_z1", 15:23))), na.rm = FALSE)
    ) / 2,
    
    ink_total_z2 = (
      rowMeans(across(all_of(sprintf("inc_%02d_z2", 1:14))), na.rm = FALSE) +
        rowMeans(across(all_of(sprintf("inc_%02d_z2", 15:23))), na.rm = FALSE)
    ) / 2,
    
    ink_total_post = (
      rowMeans(across(all_of(sprintf("inc_%02d_post", 1:14))), na.rm = FALSE) +
        rowMeans(across(all_of(sprintf("inc_%02d_post", 15:23))), na.rm = FALSE)
    ) / 2,
    
    gpd_total_pre = rowMeans(
      across(all_of(sprintf("scl_%02d_pre", 1:9))), na.rm = FALSE
    ),
    gpd_total_z1 = rowMeans(
      across(all_of(sprintf("scl_%02d_z1", 1:9))), na.rm = FALSE
    ),
    gpd_total_z2 = rowMeans(
      across(all_of(sprintf("scl_%02d_z2", 1:9))), na.rm = FALSE
    ),
    gpd_total_post = rowMeans(
      across(all_of(sprintf("scl_%02d_post", 1:9))), na.rm = FALSE
    )
  ) %>%
  filter(if_any(-patient_id, ~ !is.na(.x)))

########## Trier sample ############################################################

# INC is the mean of the ten INK-10 items.
# GPD is the mean of the eleven symptom items.
# Baseline BSI items are already transformed to the 1-4 scale.
trier_scores <- trier_items %>%
  transmute(
    patient_id,
    
    ink_total_pre = rowMeans(
      across(all_of(sprintf("inc_%02d_pre", 1:10))), na.rm = FALSE
    ),
    ink_total_s05 = rowMeans(
      across(all_of(sprintf("inc_%02d_s05", 1:10))), na.rm = FALSE
    ),
    ink_total_s10 = rowMeans(
      across(all_of(sprintf("inc_%02d_s10", 1:10))), na.rm = FALSE
    ),
    ink_total_s15 = rowMeans(
      across(all_of(sprintf("inc_%02d_s15", 1:10))), na.rm = FALSE
    ),
    ink_total_s20 = rowMeans(
      across(all_of(sprintf("inc_%02d_s20", 1:10))), na.rm = FALSE
    ),
    ink_total_s25 = rowMeans(
      across(all_of(sprintf("inc_%02d_s25", 1:10))), na.rm = FALSE
    ),
    
    gpd_total_pre = rowMeans(
      across(all_of(sprintf("gpd_%02d_pre", 1:11))), na.rm = FALSE
    ),
    gpd_total_s05 = rowMeans(
      across(all_of(sprintf("gpd_%02d_s05", 1:11))), na.rm = FALSE
    ),
    gpd_total_s10 = rowMeans(
      across(all_of(sprintf("gpd_%02d_s10", 1:11))), na.rm = FALSE
    ),
    gpd_total_s15 = rowMeans(
      across(all_of(sprintf("gpd_%02d_s15", 1:11))), na.rm = FALSE
    ),
    gpd_total_s20 = rowMeans(
      across(all_of(sprintf("gpd_%02d_s20", 1:11))), na.rm = FALSE
    ),
    gpd_total_s25 = rowMeans(
      across(all_of(sprintf("gpd_%02d_s25", 1:11))), na.rm = FALSE
    )
  ) %>%
  filter(if_any(-patient_id, ~ !is.na(.x)))

########## Analytic sample sizes ###################################################

sample_sizes <- tibble(
  sample = c("Bern", "Zurich", "Trier"),
  prepared_n = c(
    nrow(bern_items),
    nrow(zurich_items),
    nrow(trier_items)
  ),
  analytic_n = c(
    nrow(bern_scores),
    nrow(zurich_scores),
    nrow(trier_scores)
  )
) %>%
  mutate(excluded_no_model_score = prepared_n - analytic_n)

stopifnot(
  nrow(bern_scores) == 101L,
  nrow(zurich_scores) == 145L,
  nrow(trier_scores) == 1410L
)

print(sample_sizes, width = Inf)

########## Model datasets ##########################################################

# Number scores by assessment order for the model scripts.
df_bern <- bern_scores %>%
  transmute(
    ink1 = ink_total_pre,
    ink2 = ink_total_z1,
    ink3 = ink_total_z2,
    ink4 = ink_total_post,
    gpd1 = gpd_total_pre,
    gpd2 = gpd_total_z1,
    gpd3 = gpd_total_z2,
    gpd4 = gpd_total_post
  )

df_zurich <- zurich_scores %>%
  transmute(
    ink1 = ink_total_pre,
    ink2 = ink_total_z1,
    ink3 = ink_total_z2,
    ink4 = ink_total_post,
    gpd1 = gpd_total_pre,
    gpd2 = gpd_total_z1,
    gpd3 = gpd_total_z2,
    gpd4 = gpd_total_post
  )

df_trier <- trier_scores %>%
  transmute(
    ink1 = ink_total_pre,
    ink2 = ink_total_s05,
    ink3 = ink_total_s10,
    ink4 = ink_total_s15,
    ink5 = ink_total_s20,
    ink6 = ink_total_s25,
    gpd1 = gpd_total_pre,
    gpd2 = gpd_total_s05,
    gpd3 = gpd_total_s10,
    gpd4 = gpd_total_s15,
    gpd5 = gpd_total_s20,
    gpd6 = gpd_total_s25
  )

########## Export ##################################################################

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(bern_scores, "data/processed/bern_scores.csv", na = "")
write_csv(zurich_scores, "data/processed/zurich_scores.csv", na = "")
write_csv(trier_scores, "data/processed/trier_scores.csv", na = "")

save(
  df_bern,
  df_zurich,
  df_trier,
  file = "data/data.RData"
)