REPRODUCING THE SIMULATIONS AND CBASP DATA APPLICATION
====================================================

This folder contains three simulation workflows and one real-data-based
application. Each simulation script defines delta(), but the arguments and
saved outputs differ. Run each workflow in a fresh R session and keep its
results in a separate directory. The plotting reports read saved results;
they do not rerun simulations.

1. FILE GUIDE
-------------

simu_updated_empirical_bootstrap_s2_conditionalX2_gexp.R
  X-only covariate-adaptive randomization with the observed-row empirical
  bootstrap. Runs proposed and comparison procedures over three scenarios
  and four treatment effects. delta() returns a numeric length-132 vector
  of rejection indicators at its default settings.

plot_simulation_main.Rmd
  Reads the preceding script's seed files. Produces type I error/power
  summaries and proposed-versus-naive and proposed-versus-comparator figures.

simu_ball_corrected_s2_conditionalX2_gexp.R
  All-covariate balancing simulation. Proposed procedures use analytical
  variance adjustments, without bootstrap. delta() returns a length-132
  rejection-indicator vector, with a different comparator order from the
  empirical-bootstrap script.

plot_simu_ball_corrected.Rmd
  Reads ONLY the all-covariate balancing script's results, using its
  scenario-specific comparator mapping.

simu_variance_one_seed_2100_s2_conditionalX2_gexp_order3_knot0.R
  X-only empirical-bootstrap variance study over sample sizes 300, 600,
  900, 1200, 1500, 1800 and 2100. Saves adjusted treatment variances and
  independent reference treatment estimates for the two-sample and linear
  working-model procedures. delta() returns a data frame, not a length-132
  vector. The default grid produces 42 rows per seed.

variance_type1_report.Rmd
  Reads the variance-study data frames and plots treatment-variance ratios
  and empirical type I error versus sample size. Creates a reusable g2 plot.

realdata.R
  Prepares observed CBASP baseline profiles, fits observed-outcome models,
  and runs the real-data-based simulation. rd_reproduce() generates the
  four-method table and the supplementary two-sample and partially linear
  results. Executing this script with Rscript calls rd_reproduce().

CBASP.csv
  Original clinical-trial data, including A2, AGE, HAMD24 and FinalHAMD.
  Needed for realdata.R; not needed for the three synthetic simulations.

job_run.R
  Command-line worker for the empirical-bootstrap power simulation.
  Reads a seed from its first command-line argument, calls delta(seed),
  and saves object res1 as seedN.RData. Its current paths start with /home.

main_job_content.sh
  Slurm script that loads R and runs Rscript job_run.R for one seed.
  Current requests: one task, one CPU, 5G memory per CPU, 30 minutes.
  Standard output/error use the job name and job ID.

main_job.sh
  Submits one Slurm job per seed, from seed 1 through seed 1000.

Function-level purpose, arguments and return values are documented in the
R scripts. roxygen2 is not required to source these standalone scripts.

2. SOFTWARE AND WORKING DIRECTORY
--------------------------------

Use R/RStudio with the following packages. In a fresh R session, install
missing packages before loading simulation scripts:

  install.packages(
    c("Rcpp", "mgcv", "MASS", "splines2", "carat",
      "ggplot2", "patchwork", "rmarkdown", "knitr"),
    repos = "https://cloud.r-project.org"
  )

The power simulations require mgcv, MASS, splines2 and carat.
The variance worker requires mgcv, MASS and splines2; Rcpp is optional
for its faster allocation/bootstrap engine. The real-data example uses
mgcv, splines2 and carat. Reports require ggplot2, knitr and rmarkdown;
the variance report additionally requires patchwork and ggplot2 >= 3.4.0.

The carat version used with this code requires Rcpp >= 1.1.2. If an older
Rcpp is already loaded, restart R before updating it. After installation,
restart R again, then check:

  library(Rcpp)
  getNamespaceVersion("Rcpp")
  library(splines2)
  library(carat)

Do not try to unload Rcpp while splines2 imports it. If multiple R libraries
are present, check .libPaths() and find.package("Rcpp") to confirm which
installation is selected.

HTML rendering needs Pandoc. Use RStudio's Knit button, or configure Pandoc
for command-line rendering. Check rmarkdown::pandoc_available() before
running rmarkdown::render() from a standalone R installation.

Before the examples below, set the R working directory to this code folder,
for example using RStudio's Session > Set Working Directory menu. Relative
paths in source() and file.copy() refer to the working directory. During
normal knitting, report paths are resolved relative to the Rmd file.

3. CURRENT SYNTHETIC DATA-GENERATING SETTINGS
--------------------------------------------

The outcome has treatment means 1+u and 1, plus g(X)+f(Z) and zero-mean
normal noise. In these scripts u is the mean contrast for treatment 1
minus treatment 0. The nominal two-sided test level is 0.05.

Scenario (i):
  A latent bivariate normal has X-Z correlation 0.7. X is mapped into
  the discrete levels -2, -1, 1 and 2. Z is clipped to [-3,3].
  g(X)=2*X, f(Z)=exp(Z)/2.5, and the error variance is 1.

Scenario (ii):
  X is standard normal conditional on -2<X<2 and is independent of Z.
  This is continuous truncation without endpoint atoms. X is not
  standardized after sampling. Z is clipped to [-3,3].
  g(X)=exp(X), f(Z)=4*log(Z+4), and the error variance is 1.

Scenario (iii):
  The latent normal correlation matrix for (X1,X2,Z) is
    [1, 0.5, 0.7; 0.5, 1, 0; 0.7, 0, 1].
  Each X coordinate is clipped to [-1,1], and Z is clipped to [-3,3].
  g(X)=2*X1*X2^2, f(Z)=exp(Z)/2.5, and
  Var(error | X,Z)=log(1+2*X1^2+1.5*Z^2).

These are the settings currently implemented in the supplied files.
In particular, they are not the earlier Scenario-(iii) variant with
g(X)=X1*X2^2/5 and X on [-3,3]. Clipping and conditional truncation
are different distributions and their simulation results must remain
in separate batches.

The proposed CAR probability rho is 0.9. The allocation features differ
between workflows:

  Empirical-bootstrap power script: balances X. Scenario (ii) uses the
  empirical quantile-knot univariate design with default degree 3
  (order 4). Scenario (iii) uses degree 2 (order 3), zero-split local
  cells and the complete X1-by-X2 tensor spline span. The full tensor is
  formed before replacing one tensor column with an intercept, retaining
  16 features. Original features and nuisance fits are frozen during
  bootstrap; observed covariate rows are jointly resampled and CAR is
  restarted. Only sigma_m is bootstrapped.

  All-covariate script: balances X and Z. Scenarios (i)/(ii) use local
  additive order-3 splines with zero knots. Scenario (iii) uses global
  weighted balancing of the complete X1-by-X2 spline tensor and the Z
  spline basis; delta() specifies degree 3 (order 4) without interior
  knots. The proposed variance calculations do not use bootstrap.

  Variance-study script: balances X. Scenario (ii) uses FIXED degree-2
  (order-3) splines, interior knot 0 and boundaries [-2,2] for every n.
  Scenario (iii) uses the same complete order-3 tensor-span construction
  as the empirical-bootstrap design. Allocation features and nuisance
  fits remain fixed during the observed-row bootstrap.

Thus the two X-only scripts do not share the same Scenario-(ii) knot/order
settings. Keep their outputs separate even when their DGPs agree.

4. EMPIRICAL-BOOTSTRAP TYPE I ERROR AND POWER
------------------------------------------

Start a fresh R session. From the code folder, run:

  source("simu_updated_empirical_bootstrap_s2_conditionalX2_gexp.R")
  seed_dir <- "empirical_bootstrap_seeds"
  dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)

  for (seed in seq_len(1000)) {
    res1 <- delta(seed, n = 1000, b = 500, u = c(0, 0.1, 0.2, 0.3))
    save(res1, file = file.path(seed_dir, paste0("seed", seed, ".RData")))
  }

Use a fresh result directory. The full run is computationally intensive;
Section 7 describes the supplied cluster submission scripts.

At these defaults, res1 has 132 entries: 44 per scenario, in Scenario
(i), (ii), (iii) order. Entries are 0/1 rejection indicators or NA.
They are not treatment estimates or fitted variance components.
The plotting report expects all four effects in the stated order.

Copy the report beside the seed files and knit it:

  file.copy("plot_simulation_main.Rmd", seed_dir, overwrite = TRUE)
  rmarkdown::render(file.path(seed_dir, "plot_simulation_main.Rmd"))

Alternatively, open the copied report in RStudio and click Knit.
The main report has data_dir="." in its settings chunk rather than a
YAML data_dir parameter. Edit that setting if keeping the report elsewhere.

Outputs in empirical_bootstrap_seeds/figures/:
  proposed_vs_naive.png/.pdf
  proposed_vs_comparators.png/.pdf
  power_summary.csv
  method_mapping.csv
  input_manifest.csv
  missing_results_by_column.csv
The HTML report includes the type I error table at u=0.
For nonzero u, rejection rates estimate power.

5. ALL-COVARIATE BALANCING TYPE I ERROR AND POWER
----------------------------------------------

Start a new R session, then run:

  source("simu_ball_corrected_s2_conditionalX2_gexp.R")
  seed_dir <- "ball_corrected_seeds"
  dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)

  for (seed in seq_len(1000)) {
    res1 <- delta(seed, n = 1000)
    save(res1, file = file.path(seed_dir, paste0("seed", seed, ".RData")))
  }

  rmarkdown::render(
    "plot_simu_ball_corrected.Rmd",
    params = list(data_dir = normalizePath(seed_dir), object_name = "res1")
  )

The four effects are fixed inside delta(): 0, 0.1, 0.2 and 0.3.
There is no bootstrap-count argument for this delta().

Outputs in ball_corrected_seeds/figures_ball_corrected/ have the same
figure/table basenames as Section 4. The all-covariate report uses a
different comparator mapping, including a special Scenario-(ii) order.
Do not use plot_simulation_main.Rmd for these length-132 vectors.

6. TREATMENT-VARIANCE RATIOS AND TYPE I ERROR
------------------------------------------

Start a new R session, then run:

  source("simu_variance_one_seed_2100_s2_conditionalX2_gexp_order3_knot0.R")
  seed_dir <- "variance_seeds"

  for (seed in seq_len(1000)) {
    res1 <- delta(seed, b = 500,
                  n_values = seq(300, 2100, by = 300),
                  scenarios = c("i", "ii", "iii"), effect = 0)
    vc_write_rdata(res1, output_dir = seed_dir)
  }

  rmarkdown::render(
    "variance_type1_report.Rmd",
    params = list(
      results_dir = normalizePath(seed_dir),
      output_dir = "variance_type1_conditionalX2_gexp_plots"
    )
  )

vc_write_rdata() creates the result directory and writes object res1
to seedN.RData. It refuses to overwrite an existing seed by default.
Alternatively, vc_write_seed(res1, output_dir=seed_dir) writes
variance_seed_XXXX.csv. For CSV input, set the report's file_pattern to
"^variance_seed_[0-9]+[.]csv$". Choose one format per result directory
so each replicate is read only once.

delta() defaults: B=500, effect=0, base_seed=20260921,
quadrature=64, edge_grid=1001, fast=TRUE and verbose=TRUE.
The optional Rcpp engine needs a working C++ compiler. If compilation
fails, the worker falls back to R. Set fast=FALSE to request R explicitly.
The report checks backend and software/settings metadata; analyze
different backends or fitting versions as separate batches.

Each default replicate has 42 rows: three scenarios x seven sample sizes
x two methods, Proposed-t-test and Proposed-lm. Each cell saves a fitted
treatment variance and an independent reference treatment estimate.
One seed alone does not provide a Monte Carlo variance.

For each scenario/sample-size/method, the report calculates:
  variance ratio = mean(treatment_variance) / var(mc_treatment_estimate)
  test statistic = (treatment_estimate - null_value) / sqrt(treatment_variance)
  type I error = fraction of evaluable null tests rejecting at alpha=0.05.

treatment_variance already includes the factor 4/n. Do not multiply it
by 4/n again. The components base, sigma_E and sigma_m, and their sum
adjusted_sigma2, are unscaled quantities.

The ratio reference line is 1; the type I error reference line is 0.05.
Bands show pointwise Monte Carlo intervals. The numerator uncertainty
uses sd(saved variance estimates)/sqrt(number of replicates), not their
raw SD. The denominator uncertainty is the Monte Carlo SE of the
independent reference sample variance. Type I error limits are exact
binomial intervals.

Outputs in variance_type1_conditionalX2_gexp_plots/:
  variance_ratio_type1_conditionalX2_gexp.png/.pdf
  treatment_effect_variance_ratio_bands.png/.pdf
  type1_error_bands.png/.pdf
  g2.rds
  treatment_effect_variance_ratio_summary.csv
  type1_error_summary.csv
  simulation_diagnostics.csv and scenario_ii_diagnostics.csv
  treatment_effect_variance_ratio_plot_data.csv
  type1_error_plot_data.csv
  input_manifest.csv, invalid_rows.csv and scenario2_settings.csv
  analysis_settings.csv and sessionInfo.txt

Display the combined plot after knitting without rerunning simulations:

  library(ggplot2)
  library(patchwork)
  g2 <- readRDS("variance_type1_conditionalX2_gexp_plots/g2.rds")
  g2

The variance source also provides vc_summarize() for comparing fitted
bootstrap components, oracle bootstrap components and total variances
with their independent Monte Carlo references. vc_plot() exports the
corresponding level plots. These optional helpers are separate from
the variance-ratio/type-I-error figures generated by the Rmd; their
documented arguments distinguish outer rows from reference rows.

Set params$variance_y_min and params$variance_y_max to NULL together
for automatic variance axes; type1_y_max=NULL requests an automatic
type I error upper limit. These options affect display, not the summaries.
For a strictly complete analysis, set invalid_policy="stop". The default
"report" records invalid tests and gives rejection rates conditional on
evaluable tests, with counts and all-loaded bounds in the output tables.

7. RUNNING THE SUPPLIED SLURM JOBS
--------------------------------

The supplied job_run.R currently sources:
  /home/simu_updated_empirical_bootstrap_s2_conditionalX2_gexp.R
and saves results in:
  /home/simu_updated_empirical_bootstrap_s2_conditionalX2_gexp/

Before submission, edit BOTH paths in job_run.R for your cluster:
  source(): full path to the empirical-bootstrap simulation source.
  save(): writable result directory for seedN.RData.
Create the result directory before submitting jobs.

The R module must have the required packages available. Adjust the module
name and Slurm resource limits in main_job_content.sh if needed for your
cluster. Run from the directory containing both .sh files and job_run.R:

  sbatch main_job_content.sh 1

After checking the single-seed job and its log/result, submit the full batch:

  bash main_job.sh

main_job.sh submits 1000 separate jobs; it does not submit a job array.
main_job_content.sh calls Rscript job_run.R with the seed argument.
Default log filenames are new_bootstrap_<jobID>.out and
new_bootstrap_<jobID>.err in the submission directory.
These two shell scripts require Slurm; the local R examples do not.

To use the same submission wrappers for another workflow, change the
source and output paths in job_run.R. For all-covariate balancing use
delta(seed) and save its length-132 res1 vector. For the variance study,
use delta(seed) and vc_write_rdata(res1, output_dir=...), retaining the
full data frame and its metadata. Keep each workflow's batch separate.
The current job_run.R does not run realdata.R.

8. CBASP APPLICATION: TABLE AND SUPPLEMENTARY RESULTS
---------------------------------------------------

Start a new R session. Keep CBASP.csv beside realdata.R, then run:

  source("realdata.R")
  ans <- rd_reproduce()
  ans$table

This call uses 1000 replications, seeds 1:1000, noise SD 2, a fixed
order-3 age-spline knot at zero, rho=0.9 and four comparator age strata.
Or specify the settings explicitly:

  ans <- rd_reproduce(data_file = "CBASP.csv", repetitions = 1000,
                      seed_start = 1, sigma = 2, strata = 4)

CBASP treatment codes are A2=-1 for NEF, A2=0 for combination therapy,
and A2=1 for CBT alone. Excluding CBT retains 453 baseline profiles.
All 453 are used for simulated trials. Observed-outcome regressions use
the 440 complete FinalHAMD observations. Age is standardized by its
sample mean and sample SD across the 453 eligible participants.
HAMD24 is also standardized for working regressions.

Simulation outcomes use the raw baseline covariates:
  FinalHAMD = T*20.32581 + (1-T)*25.196228
              -0.972179*AGE + 0.011102*AGE^2 + 0.391615*HAMD24 + error
  error ~ N(0, 2^2).

The baseline profiles and their row order are fixed between replications;
treatment assignments and errors are newly generated for each seed.
The table reports mean treatment estimates and their empirical SD across
replications, conditional on those profiles. Its SD row is not an
average fitted standard error. No bootstrap is needed to calculate this
table's point estimates and empirical SDs.

ans$table reports Proposed-lm, Ma2015lm, Ma2022mr and Ye2018plm.
Obtain the supplementary two-sample and partially linear results from
the SAME simulation run:

  supplementary <- subset(
    ans$summary,
    method %in% c("Proposed-t-test", "Proposed-plm"),
    select = c(method, mean_estimate, empirical_sd)
  )
  print(supplementary, row.names = FALSE, digits = 4)

Use these computed values for the corresponding manuscript sentence.
Fresh results need not exactly match previously printed decimals when
seeds, allocation settings, data ordering or package versions differ.

Outputs in realdata_results/ beside realdata.R:
  Table_D1.csv: unrounded four-method table.
  Table_D1.tex: computed LaTeX table fragment; requires booktabs in the paper.
  Table_D1.html/.pdf/.png: table previews.
  all_methods_summary.csv: all six methods and Monte Carlo SEs of their means.
  supplementary_methods.csv: the two supplementary methods.
  trial_estimates.csv: all six treatment estimates for every seed.
  simulation_seeds.csv: seed sequence.
  SmoothEffect.pdf/.png: observed-outcome GAM smooth effects of AGE and HAMD24.
  observed_quadratic_coefficients.csv: actual-data calibration regression.
  observed_quadratic_fit.txt and observed_additive_fit.txt: model summaries.
  settings.txt and sessionInfo.txt: settings, input checksum and software versions.
  reproduction.rds: estimates, table, summaries, observed models and settings.

Reload the completed run instead of simulating again:

  ans <- readRDS("realdata_results/reproduction.rds")
  ans$table
  ans$summary

The table is a simulation calibrated to observed data. The observed
regression fits and SmoothEffect figure are produced from actual outcomes
and exported separately; the table is not an analysis of 1000 observed
clinical trials. Figure filenames are listed here; manuscript figure/table
numbers should be matched to the version being submitted.

Command-line alternative:

  Rscript realdata.R

This executes rd_reproduce() automatically. Available --name=value options:
  --data-file, --output-dir, --repetitions, --seed-start, --sigma,
  --strata, --results-dir.
Quote the complete argument when a path contains spaces, for example:
  Rscript realdata.R "--data-file=/path with spaces/CBASP.csv" --repetitions=1000

To read an existing compatible batch of real-data delta() results:

  source("realdata.R")
  ans <- rd_reproduce(results_dir = "/path/to/realdata_seed_files",
                      repetitions = 1000, sigma = 2)

The reader selects seed1.RData through seed1000.RData and extracts
the six reported methods from numeric length-60 vectors, normally res1.
Old vectors do not record their knot/strata settings; verify the generating
code before using them. The optional full-inference legacy workflow is:

  rd_initialize("CBASP.csv")
  res1 <- delta(1)

This real-data delta() uses noise SDs sqrt(3), 2 and sqrt(5), bootstrap
B=500 for the first two proposed methods, and four comparator strata.
It is separate from the three synthetic delta() functions above. Its
nonzero effect means rejection indicators describe power, not type I error.

9. RESULT CHECKS AND COMMON ERRORS
--------------------------------

No seed files found:
  Point the report at the directory containing completed seed files,
  not merely the source-code directory. For plot_simulation_main.Rmd,
  copy it beside the results or edit data_dir in its settings chunk.
  For the other two reports, set data_dir or results_dir as shown above.

Wrong result format:
  Empirical/all-covariate power reports need length-132 numeric vectors.
  The variance report needs variance-study data frames with trial estimates
  and variances. realdata.R's legacy reader needs length-60 vectors.
  Rejection-only files cannot reconstruct treatment variances.

Missing or invalid results:
  Inspect the reports' manifests, missing-column audits and diagnostic CSVs.
  The main power report stops on missing plotted results unless its
  allow_missing setting is changed. The all-covariate report uses available
  entries and reports counts; its allow_missing parameter is retained for
  compatibility rather than enforcing complete data. Singleton comparator
  strata or nonpositive variance estimates can produce unavailable tests.
  Such results are not counted as ordinary nonrejections.

Reproducibility:
  Keep one seed file per replication and do not mix copies across formats.
  Use a single n/effect grid, bootstrap count, algorithm and package-version
  batch for each report. The two rejection-vector formats cannot identify
  their generating settings automatically; retain the source and run settings.
  The variance study and real-data example save explicit settings/seed records.
  Power reports read all matching seed files, not a built-in seed-1:2000 subset.
  Their file manifests show which files were actually included.

Figures without rerunning simulations:
  Knit the matching Rmd using existing seed results. For the variance
  figure, reload g2.rds. For the real-data table/results, reload
  reproduction.rds. Raw results remain unchanged by plotting.
