Data documentation
==================

Scope
-----

This document describes the prepared input datasets, analytic score files,
and model objects used for the Bern, Zurich, and Trier samples.

The documentation refers to:

- scripts/preparation/01_build_analysis_data.R
- scripts/descriptives/01_sample_descriptives.R
- scripts/descriptives/02_wave_descriptives.R
- scripts/descriptives/03_baseline_reliability.R


Data organisation
-----------------

All datasets use wide format, with one row per patient. patient_id is
sample-specific.

Inputs

File                                    Rows   Contents
-------------------------------------  -----   ----------------------------
data_raw/bern_master_inc_gpd.csv          104   Participant characteristics,
                                               diagnoses, INK-23 and
                                               BSI-53 items

data_raw/zurich_master_inc_gpd.csv        149   Participant characteristics,
                                               diagnoses, INK-23 and
                                               SCL-K-9 items

data_raw/trier_master_inc_gpd.csv       1,410   Participant characteristics,
                                               diagnoses, INK-10 and
                                               symptom items

The input datasets contain prepared item responses. Required item
reversals and the Trier baseline symptom-item transformation have
already been applied.

The Trier input excludes patients with a documented age below 18.

Generated files

01_build_analysis_data.R creates:

- data/processed/bern_scores.csv
- data/processed/zurich_scores.csv
- data/processed/trier_scores.csv
- data/data.RData

The score CSVs retain patient_id and the total scores.
data.RData contains only the model-score objects.


Participant characteristics
---------------------------

Variable        Samples        Definition
--------------  -------------  --------------------------------------------
patient_id      All            Sample-specific patient identifier
age_patient     All            Age in years
sex_patient     All            Recorded category, coded as shown below
treatment_arm   Bern, Zurich   Trial treatment assignment
n_sessions      All            Total number of received treatment sessions;
                               excludes booster sessions in Bern

Recorded sex categories

Sample   Coding
-------  ---------------------------------------------------------------
Bern     M = male; W = female
Zurich   m = male; w = female
Trier    1 = male; 2 = female; 3 = diverse

Treatment arms

- Bern: EFT, SR.
- Zurich: CBT, ECBT.
- Trier: no treatment-arm variable is used.


Diagnosis variables
-------------------

Sample   Variable                       Definition
-------  -----------------------------  -----------------------------------
Bern     diagnosis_dsm_pre              Primary diagnosis in DSM-IV coding
Bern     diagnosis_icd_pre              Primary diagnosis in ICD coding
Bern     diagnosis_additional_dsm_pre    First additional diagnosis in
                                        DSM-IV coding
Bern     diagnosis_additional_icd_pre    First additional diagnosis in
                                        ICD coding
Zurich   diagnosis_dsm_pre              Primary-diagnosis text; MDE denotes
                                        a major depressive episode
Zurich   comorbid_axis1                 Additional DSM-IV Axis I diagnosis:
                                        0 = no, 1 = yes, NA = missing
Zurich   comorbid_axis2                 Additional DSM-IV Axis II diagnosis:
                                        0 = no, 1 = yes, NA = missing
Trier    diagnosis_1                    Primary diagnosis recorded using
                                        ICD-10 codes
Trier    diagnosis_2 to diagnosis_5     Additional diagnoses recorded using
                                        ICD-10 codes

The descriptive script standardises whitespace and treats empty
diagnosis strings as missing.

Primary-diagnosis grouping

Bern

Primary diagnoses are grouped using the following recorded DSM-IV codes:

- Depressive disorders:
  296.2X, 296.3X, 300.4, 311.

- Anxiety or adjustment disorders:
  300.00, 300.01, 300.02, 300.21, 300.23, 300.29, 309.81,
  309.0, 309.24, 309.28, 309.9.

The script extracts the first code and matches it exactly, ignoring
letter case.

Zurich

Primary-diagnosis entries containing MDE as a separate word are classified
as a documented major depressive episode.

Trier

diagnosis_1 is grouped by its ICD-10 F-code range:

Range     Group
--------  ---------------------------------------------------------------
F00-F09   Organic mental disorders
F10-F19   Substance-related disorders
F20-F29   Schizophrenia, schizotypal and delusional disorders
F30-F39   Mood/affective disorders
F40-F48   Neurotic, stress-related and somatoform disorders
F50-F59   Behavioural syndromes associated with physiological disturbances
F60-F69   Disorders of adult personality and behaviour
F70-F79   Intellectual disability
F80-F89   Disorders of psychological development
F90-F98   Behavioural and emotional disorders with onset in childhood
F99       Unspecified mental disorder

For all samples, the script output distinguishes missing primary
diagnoses from entries that cannot be assigned using these rules.

Additional diagnoses

For Bern, an additional diagnosis is recorded as present if either
additional-diagnosis field is nonmissing.

For Zurich:

- Yes: at least one axis equals 1, regardless of the other axis.
- No: both axes equal 0.
- Unknown: neither axis equals 1, and at least one axis is missing.

For Trier, an additional diagnosis is recorded as present if any of
diagnosis_2 to diagnosis_5 is nonmissing.

For Bern and Trier, empty additional-diagnosis fields are treated as
no documented additional diagnosis.

All diagnosis percentages use the full analytic sample as denominator,
including patients with missing or unclear primary diagnoses.


Item variables and assessments
------------------------------

Item names contain a prefix, a two-digit item number, and an assessment
suffix. For example, inc_01_pre is the first prepared INC item at
pretreatment.

Sample   Item prefixes and numbers
-------  ---------------------------------------------------------------
Bern     inc_01 to inc_23; bsi_01 to bsi_53
Zurich   inc_01 to inc_23; scl_01 to scl_09
Trier    inc_01 to inc_10; gpd_01 to gpd_11

Sample   Suffix   Assessment
-------  -------  ------------------------------------------------------
Bern     pre      Pretreatment
         z1       Session 8
         z2       Session 16
         post     Posttreatment

Zurich   pre      Pretreatment
         z1       Session 7
         z2       Session 14
         post     Posttreatment

Trier    pre      Pretreatment
         s05      Session 5
         s10      Session 10
         s15      Session 15
         s20      Session 20
         s25      Session 25


Total-score calculation
-----------------------

Every score requires complete responses to all contributing items.
A missing constituent item produces a missing total. Other scores for
the same patient are calculated independently.

INC scoring

All INC items are supplied in the scoring direction used for the analyses.
No item reversal is performed by the build script.

Bern and Zurich

The INK-23 total gives equal weight to the approach and avoidance
subscale means:

    INC = (mean(items 1-14) + mean(items 15-23)) / 2

It is not the simple mean of all 23 items.

Trier

INC is the mean of all ten prepared INK-10 items.

INC scores range from 1 to 5 in all samples. Higher values indicate
greater incongruence.

GPD: Bern and Zurich

- Bern: mean of all 53 BSI items.
- Zurich: mean of all nine SCL-K-9 items.

GPD: Trier

GPD is the mean of eleven prepared symptom items.

At pretreatment, the distress score was calculated from the following 11 BSI variables: BSI001, BSI016, BSI017, BSI018, BSI019, BSI025, BSI035, BSI038, BSI039, BSI045, and BSI050.

The baseline values have already been transformed using y = 1 + 0.75x, 
yielding a possible range of 1–4. No further transformation is performed.

At sessions 5-25, the eleven items are taken from HSC001 to HSC011.

Possible GPD scores range from 0 to 4 in Bern and Zurich and from 1-4
at all assessments in Trier.

Higher scores indicate greater psychological distress.


Analytic sample selection
-------------------------

The build script retains patients with at least one observed INC or GPD
total across the modelled assessments.

Sample   Prepared N   Analytic N   Excluded: all model scores missing
-------  -----------  -----------  ----------------------------------
Bern             104          101                                  3
Zurich           149          145                                  4
Trier          1,410        1,410                                  0

01_sample_descriptives.R and 03_baseline_reliability.R use the saved
score-file IDs to select the same analytic patients from the prepared
item datasets. Missing demographic or diagnostic information does not
cause further exclusions.


Analytic score files and model objects
-------------------------------------

The score CSVs contain patient_id, ink_total_<assessment>, and
gpd_total_<assessment>. Numerical values are not rounded before export.

data/data.RData contains:

Object       Rows   Score variables
-----------  -----  ----------------------------------------------------
df_bern        101  ink1 to ink4, gpd1 to gpd4
df_zurich      145  ink1 to ink4, gpd1 to gpd4
df_trier     1,410  ink1 to ink6, gpd1 to gpd6

The numbers indicate the chronological assessment order listed above.
These model objects contain no patient IDs or participant characteristics.


Descriptive statistics and reliability
--------------------------------------

Age and session-count means and SDs use available observations for
each variable. Categorical percentages use the full analytic sample.
The script output includes available and missing counts.

Wave means and SDs use available total scores separately for each
construct and assessment; complete data across all waves are not required.

Raw Cronbach's alpha uses complete baseline item sets, selected
separately for each instrument within the analytic sample:

- Bern: INK-23 and BSI-53.
- Zurich: INK-23 and SCL-K-9.
- Trier: INK-10 and the eleven mapped baseline BSI items.

Automatic item reversal and removal are disabled.

For INK-23, alpha describes ordinary item-level internal consistency.
It is not a separate reliability estimate tailored to the equal-weighted
combination of the two subscale means.

Rounding is applied only to displayed results.