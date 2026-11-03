Project overview

This project estimates random intercept cross-lagged panel models
(RI-CLPMs) between motivational incongruence (INC) and general
psychological distress (GPD), for the Bern, Zurich, and Trier samples.

Because participants went through an intervention, the sample-level
means of INC and GPD are expected to decrease across measurement 
waves as the intervention progresses. Whether "to detrend or not to detrend" 
such data (and how) is addressed by Falkenström et al. (2026). Since we 
use the RI-CLPM (Hamaker et al., 2015), we follow the recommendations 
in both papers: first to address detrending, and second to estimate the 
prospective effects between INC and GPD in all three samples.

The workflow begins with three separately prepared item-level datasets
for Bern, Zurich, and Trier. The build script calculates the analysis
scores; three descriptive scripts calculate the manuscript's descriptive
statistics. These scripts are described below, followed by the models.
Participant-level data are not publicly shared. Running the scripts
requires access to the prepared datasets.

The analyses follow 5 steps. Each has its own script and output folders;
each answers the question the previous step raises:

  1. Mean structure (01 mean structure/): is there a structural change
     (a trend) in the wave means, tested jointly for INC and GPD?
  2. RI-CLPM mean structure (02 riclpm mean structure/): given a trend
     exists, does an RI-CLPM whose intercepts are free to differ across
     waves ("detrended") fit better than one that assumes a single,
     constant mean ("nondetrended"; analogue to Falkenström et al., 2026)?
  3. RI-CLPM (03 riclpm/): using the detrended RI-CLPM from step 2, can the
     autoregressive and cross-lagged effects be constrained to be equal across
     waves (Figure 1)?
  4. LCM-SR exploratory analysis (04 lcm-sr exploratory analysis/): a
     separate check on step 3's conclusions using a latent curve model
     with structured residuals (Curran et al., 2014), which models
     between-person differences in the rate of change.
  5. Mean trend figure (05_mean_trend_figure/): combines step 1's three
     per-sample plots into one figure.

Run all scripts with the project root as the working directory. First run:

  source("scripts/preparation/01_build_analysis_data.R")

Then run the three scripts in scripts/descriptives/ and the model scripts
in scripts/models/. Each descriptive script reads its own inputs and
prints its results to the R console. All five model steps use the same
prepared scores from data/data.RData, directly or through a sourced script.
Steps 1, 3, and 4 are self-contained and can be run directly, for example:

  source("scripts/models/01_mean_structure.R")

Step 2 is the exception: its "detrended" model is exactly step 3's RI-CLPM,
so rather than fitting it a second time, 02_riclpm_mean_structure.R
sources scripts/models/03_riclpm.R and reuses its fitted models. Running
step 2 therefore also (re)produces step 3's output; this is intentional,
so its detrended results and plots remain in the step-3 folders. Step 5
sources scripts/models/01_mean_structure.R and regenerates step 1
before producing the combined figure.

Each of steps 1-4 runs in three parts: fit the models, save each model's
raw R/lavaan output as a text file in results/, and save plots as PNGs in
plots/. All in the folder matching the script's own number. Step 5 adds
only a plot; its sourced step 1 also writes results.

Folders

scripts/preparation/

- 01_build_analysis_data.R
  Reads the three prepared item files from data_raw/. Required item
  reversals and transformations have already been applied.
  INC gives equal weight to the approach and avoidance subscale means in
  Bern and Zurich, and is the mean of ten items in Trier. GPD is the mean
  of 53 BSI items in Bern, nine SCL-K-9 items in Zurich, and eleven prepared
  symptom items in Trier. Each score requires complete contributing items.
  Patients with at least one observed score across the modelled assessments
  are retained. The script saves one score CSV per sample with patient IDs
  in data/processed/, and the score-only objects df_bern, df_zurich, and
  df_trier in data/data.RData.

scripts/descriptives/

- 01_sample_descriptives.R
  Calculates age, session counts, sex-category and treatment-arm
  frequencies, primary-diagnosis groups, and documented additional
  diagnoses. Uses saved analytic IDs to select patients from the prepared
  item files. Continuous statistics use available observations; categorical
  percentages use the full analytic sample, with missing information
  reported separately.

- 02_wave_descriptives.R
  Reads the saved score CSVs and calculates available-case means, SDs,
  and sample sizes for INC and GPD at each assessment, producing the
  manuscript's descriptive table.

- 03_baseline_reliability.R
  Calculates raw Cronbach's alpha separately for each instrument and
  sample, using complete baseline item sets within the saved analytic
  sample. No items are automatically reversed or removed.

scripts/models/
Contains the five R scripts described above, one per step.

- 01_mean_structure.R
  For each sample, jointly tests whether INC means are equal across waves
  and GPD means are equal across waves. Every variance and covariance
  among the wave scores is free (a saturated covariance structure).
  A significant joint test indicates that at least one construct's means
  differ; it does not establish change in both constructs separately.
  Builds each sample's FIML-estimated mean plot with 95% confidence bars
  (combined into one figure in step 5). These estimates can differ from
  the available-case means in 02_wave_descriptives.R when data are missing.

- 02_riclpm_mean_structure.R
  Sources 03_riclpm.R for its detrended RI-CLPM (each wave's intercept
  free; the model used in step 3), and fits a new nondetrended RI-CLPM
  (each construct's intercept held equal across waves) alongside
  it, with the usual unconstrained/constrained (equal-over-time
  autoregressive/cross-lagged paths) pair and RMSEA_D index.

- 03_riclpm.R
  The project's main RI-CLPM (Bern/Zurich: 4 waves; Trier: 6 waves), with
  the random-intercept loading fixed to 1 and each wave's own intercept
  estimated freely (the detrended specification compared in step 2). Fits
  unconstrained and constrained (equal-over-time) autoregressive/
  cross-lagged paths and tests between them with RMSEA_D. The constrained
  version is the model in Figure 1. Summary files also include standardized
  paths and the random-intercept correlation with 95% confidence intervals.

- 04_lcmsr_exploratory_analysis.R
  Latent curve model with structured residuals (LCM-SR, Curran et al.,
  2014): each construct gets a full latent growth curve (intercept and
  slope), so between-person differences in the rate of change are captured.
  The autoregressive and cross-lagged paths then run among what is left
  once that growth trend is removed (Falkenström et al., 2026). Kept
  separate as an exploratory check on step 3's conclusions, since it makes
  a stronger demand on the data (see the admissibility check in its
  output).

- 05_mean_trend_figure.R
  Sources step 1 and stacks its three per-sample FIML-estimated mean plots
  (Bern, Zurich, Trier) into one combined figure.

data_raw/
Required local prepared item datasets: bern_master_inc_gpd.csv,
zurich_master_inc_gpd.csv, and trier_master_inc_gpd.csv.

data/
Contains generated analysis data and the data documentation.

- README.txt
  Input variables, scoring, missing-data rules, and generated datasets.
- processed/
  Generated bern_scores.csv, zurich_scores.csv, and trier_scores.csv,
  retaining patient IDs and total scores.

- data.RData
  Bern, Zurich, and Trier samples: INC and GPD scale scores per wave.

plots/
Each numbered subfolder holds that step's plots: 01 mean structure/ has
one FIML-estimated mean plot per sample; 02 riclpm mean structure/ has the
nondetrended RI-CLPM's path diagrams; 03 riclpm/ has the main analysis's
path diagrams (Figure 1, also the "detrended" side of step 2's
comparison); 04 lcm-sr exploratory analysis/ has the LCM-SR's path
diagrams; 05_mean_trend_figure/ has the combined, three-sample version of
step 1's plots.

results/
Text output straight from each fitted model, exactly as R and lavaan
produce it. Steps 2-4 split their output into two subfolders:

- model summaries/ --> one file per model (unconstrained or constrained):
  its full parameter estimates and fit indices; step 3 also includes
  standardized effects with 95% confidence intervals.
- model comparisons/ --> one file per comparison between two models: the
  chi-square difference test and RMSEA_D. Step 4 also has
  lcm_sr_admissibility_check.txt here (which LCM-SR models did and did not
  reach an admissible solution), since it compares admissibility across
  models rather than summarizing a single one.

Step 1 is not split this way: each sample's single file already combines
a model summary with its own built-in test (see above), so there is
nothing separate to compare it against. Step 5 has no results/ output --
it adds only a plot, while sourcing step 1 regenerates step-1 results.

Directory structure

├── README.txt
├── data
│   ├── README.txt
│   ├── data.RData
│   └── processed
│       ├── bern_scores.csv
│       ├── trier_scores.csv
│       └── zurich_scores.csv
├── data_raw
│   ├── bern_master_inc_gpd.csv
│   ├── trier_master_inc_gpd.csv
│   └── zurich_master_inc_gpd.csv
├── plots
│   ├── 01 mean structure
│   │   ├── mean_trend_bern.png
│   │   ├── mean_trend_trier.png
│   │   └── mean_trend_zurich.png
│   ├── 02 riclpm mean structure
│   │   ├── riclpm_nondetrended_bern_est_std.png
│   │   ├── riclpm_nondetrended_bern_names.png
│   │   ├── riclpm_nondetrended_trier_est_std.png
│   │   ├── riclpm_nondetrended_trier_names.png
│   │   ├── riclpm_nondetrended_zurich_est_std.png
│   │   └── riclpm_nondetrended_zurich_names.png
│   ├── 03 riclpm
│   │   ├── riclpm_bern_est_std.png
│   │   ├── riclpm_bern_names.png
│   │   ├── riclpm_trier_est_std.png
│   │   ├── riclpm_trier_names.png
│   │   ├── riclpm_zurich_est_std.png
│   │   └── riclpm_zurich_names.png
│   ├── 04 lcm-sr exploratory analysis
│   │   ├── lcm_sr_bern_est_std.png
│   │   ├── lcm_sr_bern_names.png
│   │   ├── lcm_sr_trier_est_std.png
│   │   ├── lcm_sr_trier_names.png
│   │   ├── lcm_sr_zurich_est_std.png
│   │   └── lcm_sr_zurich_names.png
│   └── 05_mean_trend_figure
│       └── mean_trend_combined_new.png
├── results
│   ├── 01 mean structure
│   │   ├── bern_4w_mean_structure.txt
│   │   ├── trier_6w_mean_structure.txt
│   │   └── zurich_4w_mean_structure.txt
│   ├── 02 riclpm mean structure
│   │   ├── model comparisons
│   │   │   ├── bern_4w_detrending_comparison.txt
│   │   │   ├── nondetrended_bern_4w_model_comparison.txt
│   │   │   ├── nondetrended_trier_6w_model_comparison.txt
│   │   │   ├── nondetrended_zurich_4w_model_comparison.txt
│   │   │   ├── trier_6w_detrending_comparison.txt
│   │   │   └── zurich_4w_detrending_comparison.txt
│   │   └── model summaries
│   │       ├── nondetrended_bern_4w_constrained_summary.txt
│   │       ├── nondetrended_bern_4w_unconstrained_summary.txt
│   │       ├── nondetrended_trier_6w_constrained_summary.txt
│   │       ├── nondetrended_trier_6w_unconstrained_summary.txt
│   │       ├── nondetrended_zurich_4w_constrained_summary.txt
│   │       └── nondetrended_zurich_4w_unconstrained_summary.txt
│   ├── 03 riclpm
│   │   ├── model comparisons
│   │   │   ├── bern_4w_model_comparison.txt
│   │   │   ├── trier_6w_model_comparison.txt
│   │   │   └── zurich_4w_model_comparison.txt
│   │   └── model summaries
│   │       ├── bern_4w_constrained_summary.txt
│   │       ├── bern_4w_unconstrained_summary.txt
│   │       ├── trier_6w_constrained_summary.txt
│   │       ├── trier_6w_unconstrained_summary.txt
│   │       ├── zurich_4w_constrained_summary.txt
│   │       └── zurich_4w_unconstrained_summary.txt
│   └── 04 lcm-sr exploratory analysis
│       ├── model comparisons
│       │   ├── lcm_sr_admissibility_check.txt
│       │   ├── lcm_sr_bern_4w_model_comparison.txt
│       │   ├── lcm_sr_trier_6w_model_comparison.txt
│       │   └── lcm_sr_zurich_4w_model_comparison.txt
│       └── model summaries
│           ├── lcm_sr_bern_4w_constrained_summary.txt
│           ├── lcm_sr_bern_4w_unconstrained_summary.txt
│           ├── lcm_sr_trier_6w_constrained_summary.txt
│           ├── lcm_sr_trier_6w_unconstrained_summary.txt
│           ├── lcm_sr_zurich_4w_constrained_summary.txt
│           └── lcm_sr_zurich_4w_unconstrained_summary.txt
└── scripts
    ├── descriptives
    │   ├── 01_sample_descriptives.R
    │   ├── 02_wave_descriptives.R
    │   └── 03_baseline_reliability.R
    ├── models
    │   ├── 01_mean_structure.R
    │   ├── 02_riclpm_mean_structure.R
    │   ├── 03_riclpm.R
    │   ├── 04_lcmsr_exploratory_analysis.R
    │   └── 05_mean_trend_figure.R
    └── preparation
        └── 01_build_analysis_data.R

References

Curran, P. J., Howard, A. L., Bainter, S. A., Lane, S. T., & McGinley, J. S.
(2014). The separation of between-person and within-person components of
individual change over time: A latent curve model with structured
residuals. Journal of Consulting and Clinical Psychology, 82(5), 879-894.
https://doi.org/10.1037/a0035297

Falkenström, F., Solomonov, N., & Rubel, J. (2026). To detrend, or not to
detrend, that is the question? The effects of detrending on cross-lagged
effects in panel models. Psychological Methods, 31(2), 252–270.
https://doi.org/10.1037/met0000632

Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique of
the cross-lagged panel model. Psychological Methods, 20(1), 102-116.
https://doi.org/10.1037/a0038889
