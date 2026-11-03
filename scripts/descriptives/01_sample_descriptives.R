# Calculates participant characteristics and diagnosis distributions.
# Percentages use the full analytic sample as the denominator.
# Run from the project root directory after building the analysis datasets.

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

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

# Use the saved analytic IDs without excluding missing demographic information.
bern_analysis <- bern_scores %>%
  select(patient_id) %>%
  left_join(bern_items, by = "patient_id")

zurich_analysis <- zurich_scores %>%
  select(patient_id) %>%
  left_join(zurich_items, by = "patient_id")

trier_analysis <- trier_scores %>%
  select(patient_id) %>%
  left_join(trier_items, by = "patient_id")

########## Participant characteristics #############################################

# Retain sex categories as coded in the prepared source datasets.
participant_characteristics <- bind_rows(
  bern_analysis %>%
    transmute(
      sample = "Bern",
      age = age_patient,
      sessions = n_sessions,
      sex_source = as.character(sex_patient),
      treatment_arm = as.character(treatment_arm)
    ),
  
  zurich_analysis %>%
    transmute(
      sample = "Zurich",
      age = age_patient,
      sessions = n_sessions,
      sex_source = as.character(sex_patient),
      treatment_arm = as.character(treatment_arm)
    ),
  
  trier_analysis %>%
    transmute(
      sample = "Trier",
      age = age_patient,
      sessions = n_sessions,
      sex_source = as.character(sex_patient),
      treatment_arm = NA_character_
    )
) %>%
  mutate(across(
    c(sex_source, treatment_arm),
    ~ na_if(str_squish(.x), "")
  ))

########## Age and sessions ########################################################

# Means and SDs use available observations for each variable.
continuous_summary <- participant_characteristics %>%
  select(sample, age, sessions) %>%
  pivot_longer(
    cols = c(age, sessions),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(sample, variable) %>%
  summarise(
    analytic_n = n(),
    available_n = sum(!is.na(value)),
    missing_n = sum(is.na(value)),
    mean = if (any(!is.na(value))) {
      mean(value, na.rm = TRUE)
    } else {
      NA_real_
    },
    sd = if (sum(!is.na(value)) > 1L) {
      sd(value, na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

########## Sex categories and treatment arms #######################################

categorical_values <- participant_characteristics %>%
  select(sample, sex_source, treatment_arm) %>%
  pivot_longer(
    cols = c(sex_source, treatment_arm),
    names_to = "variable",
    values_to = "category"
  ) %>%
  filter(!(sample == "Trier" & variable == "treatment_arm"))

categorical_availability <- categorical_values %>%
  group_by(sample, variable) %>%
  summarise(
    analytic_n = n(),
    available_n = sum(!is.na(category)),
    missing_n = sum(is.na(category)),
    .groups = "drop"
  )

categorical_summary <- categorical_values %>%
  count(sample, variable, category, name = "n") %>%
  left_join(
    categorical_availability,
    by = c("sample", "variable")
  ) %>%
  mutate(
    category = coalesce(category, "Missing"),
    percent_analytic = 100 * n / analytic_n
  ) %>%
  select(
    sample, variable, category,
    analytic_n, available_n, missing_n, n, percent_analytic
  )

########## Diagnosis fields ########################################################

# Empty diagnosis fields are treated as missing.
bern_analysis <- bern_analysis %>%
  mutate(across(
    c(
      diagnosis_dsm_pre,
      diagnosis_additional_dsm_pre,
      diagnosis_additional_icd_pre
    ),
    ~ na_if(str_squish(as.character(.x)), "")
  ))

zurich_analysis <- zurich_analysis %>%
  mutate(across(
    c(diagnosis_dsm_pre, comorbid_axis1, comorbid_axis2),
    ~ na_if(str_squish(as.character(.x)), "")
  ))

trier_analysis <- trier_analysis %>%
  mutate(across(
    all_of(paste0("diagnosis_", 1:5)),
    ~ na_if(str_squish(as.character(.x)), "")
  ))

stopifnot(
  all(
    is.na(zurich_analysis$comorbid_axis1) |
      zurich_analysis$comorbid_axis1 %in% c("0", "1")
  ),
  all(
    is.na(zurich_analysis$comorbid_axis2) |
      zurich_analysis$comorbid_axis2 %in% c("0", "1")
  )
)

########## Bern primary diagnoses ##################################################

bern_diagnoses <- bern_analysis %>%
  transmute(
    sample = "Bern",
    diagnosis_original = diagnosis_dsm_pre,
    diagnosis_code = str_to_upper(word(diagnosis_dsm_pre, 1))
  ) %>%
  mutate(
    diagnosis_group = case_when(
      is.na(diagnosis_original) ~ "Missing diagnosis",
      
      diagnosis_code %in% c("296.2X", "296.3X", "300.4", "311") ~
        "Depressive disorders",
      
      diagnosis_code %in% c(
        "300.00", "300.01", "300.02",
        "300.21", "300.23", "300.29", "309.81",
        "309.0", "309.24", "309.28", "309.9"
      ) ~ "Anxiety- or adjustment disorders",
      
      TRUE ~ "Unclear / unmapped entry"
    )
  )

########## Zurich primary diagnoses ################################################

zurich_diagnoses <- zurich_analysis %>%
  transmute(
    sample = "Zurich",
    diagnosis_original = diagnosis_dsm_pre
  ) %>%
  mutate(
    diagnosis_group = case_when(
      is.na(diagnosis_original) ~ "Missing diagnosis",
      str_detect(str_to_upper(diagnosis_original), "\\bMDE\\b") ~
        "Documented major depressive episode",
      TRUE ~ "Unclear / unmapped entry"
    )
  )

########## Trier primary diagnoses #################################################

# Use diagnosis_1 as the primary diagnosis.
# Assign broad ICD-10 groups only to entries with an F-code structure.
trier_diagnoses <- trier_analysis %>%
  transmute(
    sample = "Trier",
    diagnosis_original = diagnosis_1,
    diagnosis_code = str_to_upper(diagnosis_1)
  ) %>%
  mutate(
    code_structure_ok = coalesce(
      str_detect(diagnosis_code, "^F[0-9]{2}(\\.[0-9]+)?$"),
      FALSE
    ),
    icd_number = as.integer(if_else(
      code_structure_ok,
      str_sub(diagnosis_code, 2, 3),
      NA_character_
    )),
    diagnosis_group = case_when(
      is.na(diagnosis_original) ~ "Missing diagnosis",
      !code_structure_ok ~ "Unclear / unmapped entry",
      
      between(icd_number, 0, 9) ~
        "F00-F09 Organic mental disorders",
      between(icd_number, 10, 19) ~
        "F10-F19 Substance-related disorders",
      between(icd_number, 20, 29) ~
        "F20-F29 Schizophrenia, schizotypal and delusional disorders",
      between(icd_number, 30, 39) ~
        "F30-F39 Mood / affective disorders",
      between(icd_number, 40, 48) ~
        "F40-F48 Neurotic, stress-related and somatoform disorders",
      between(icd_number, 50, 59) ~
        "F50-F59 Behavioural syndromes associated with physiological disturbances",
      between(icd_number, 60, 69) ~
        "F60-F69 Disorders of adult personality and behaviour",
      between(icd_number, 70, 79) ~
        "F70-F79 Intellectual disability",
      between(icd_number, 80, 89) ~
        "F80-F89 Disorders of psychological development",
      between(icd_number, 90, 98) ~
        "F90-F98 Behavioural and emotional disorders with onset in childhood",
      icd_number == 99 ~
        "F99 Unspecified mental disorder",
      
      TRUE ~ "Unclear / unmapped entry"
    )
  )

########## Primary-diagnosis summary ###############################################

primary_diagnoses <- bind_rows(
  bern_diagnoses %>%
    select(sample, diagnosis_original, diagnosis_group),
  zurich_diagnoses %>%
    select(sample, diagnosis_original, diagnosis_group),
  trier_diagnoses %>%
    select(sample, diagnosis_original, diagnosis_group)
)

diagnosis_availability <- primary_diagnoses %>%
  group_by(sample) %>%
  summarise(
    analytic_n = n(),
    recorded_n = sum(!is.na(diagnosis_original)),
    missing_n = sum(is.na(diagnosis_original)),
    unclear_n = sum(diagnosis_group == "Unclear / unmapped entry"),
    classifiable_n = sum(
      !diagnosis_group %in% c(
        "Missing diagnosis",
        "Unclear / unmapped entry"
      )
    ),
    .groups = "drop"
  )

primary_summary <- primary_diagnoses %>%
  count(sample, diagnosis_group, name = "n") %>%
  left_join(diagnosis_availability, by = "sample") %>%
  mutate(percent_analytic = 100 * n / analytic_n) %>%
  select(
    sample, diagnosis_group, analytic_n,
    recorded_n, missing_n, unclear_n, classifiable_n,
    n, percent_analytic
  )

########## Additional diagnoses ####################################################

# Bern: either DSM or ICD can document the first additional diagnosis.
# Empty fields indicate no documented additional diagnosis.
bern_additional <- bern_analysis %>%
  transmute(
    sample = "Bern",
    category = if_else(
      !is.na(diagnosis_additional_dsm_pre) |
        !is.na(diagnosis_additional_icd_pre),
      "Yes",
      "No"
    )
  )

# Zurich: either positive axis establishes an additional diagnosis.
# Both axes must be negative to classify the status as No.
zurich_additional <- zurich_analysis %>%
  transmute(
    sample = "Zurich",
    category = case_when(
      comorbid_axis1 == "1" | comorbid_axis2 == "1" ~ "Yes",
      comorbid_axis1 == "0" & comorbid_axis2 == "0" ~ "No",
      TRUE ~ "Unknown"
    )
  )

# Trier: count any recorded entry in the additional diagnosis fields.
# This indicator does not count distinct diagnoses.
trier_additional <- trier_analysis %>%
  transmute(
    sample = "Trier",
    category = if_else(
      if_any(all_of(paste0("diagnosis_", 2:5)), ~ !is.na(.x)),
      "Yes",
      "No"
    )
  )

additional_values <- bind_rows(
  bern_additional,
  zurich_additional,
  trier_additional
)

additional_availability <- additional_values %>%
  group_by(sample) %>%
  summarise(
    analytic_n = n(),
    known_n = sum(category != "Unknown"),
    unknown_n = sum(category == "Unknown"),
    .groups = "drop"
  )

additional_summary <- additional_values %>%
  count(sample, category, name = "n") %>%
  group_by(sample) %>%
  complete(
    category = c("Yes", "No", "Unknown"),
    fill = list(n = 0L)
  ) %>%
  ungroup() %>%
  left_join(additional_availability, by = "sample") %>%
  mutate(percent_analytic = 100 * n / analytic_n) %>%
  select(
    sample, category, analytic_n,
    known_n, unknown_n, n, percent_analytic
  )

########## Output ##################################################################

# Format displayed values only.
print(
  continuous_summary %>%
    mutate(across(c(mean, sd), ~ sprintf("%.2f", .x))),
  n = Inf,
  width = Inf
)

print(
  categorical_summary %>%
    mutate(percent_analytic = sprintf("%.2f", percent_analytic)),
  n = Inf,
  width = Inf
)

print(
  primary_summary %>%
    mutate(percent_analytic = sprintf("%.2f", percent_analytic)),
  n = Inf,
  width = Inf
)

print(
  additional_summary %>%
    mutate(percent_analytic = sprintf("%.2f", percent_analytic)),
  n = Inf,
  width = Inf
)