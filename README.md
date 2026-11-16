# Sharp minimax rates and adaptive estimation of latent proportions in bivariate Poisson mixtures

This repository provides the R scripts, simulation procedures, and archived numerical results accompanying the manuscript entitled Sharp minimax rates and adaptive estimation of latent proportions in bivariate Poisson mixtures by Essoham Ali.

## Repository structure

- `scripts/ejs_simulations.R`: main simulation engine (base R).
- `scripts/figure4_article.R`: **manuscript Figure 4** (practical ratio-selection and paired relative MSE with 2,000 bootstrap resamples). It takes the **CSV file path** followed by the **output directory**.
- `scripts/table3_degree_sensitivity.R`: additional degree sensitivity experiment.
- `scripts/table4_support_sensitivity.R`: independent support sensitivity audit.
- `results/reference/`: archived Monte Carlo results (including the full 1,000-replication data), reference figures 1–5, including the updated manuscript Figure 4 (`figure4.pdf`).
- `results/reference/local_full/`: archived local-alternative experiment.
- `results/generated/table4_support_sensitivity/full/`: retained independent Table 4 audit.

The old two-panel theoretical/practical selection graphic formerly named `figure4_adaptive_ratio_selection` is intentionally **excluded**: it is not manuscript Figure 4. Old scripts, pilot runs, duplicate `generated/full` data, and obsolete graphic variants were also excluded.

## Reproduction commands (from repository root)

All scripts use **base R**; no Python or external R packages are required. Run the commands from the repository root, with `Rscript` available on your PATH. In RStudio, use `system2("Rscript", ...)` instead of `source()` for the CLI driver, because its test mode calls `quit()`.

```sh
# Deterministic unit checks (fast; no Monte Carlo study)
Rscript scripts/ejs_simulations.R test

# Main simulation: replication data, summaries, and manuscript Figures 1–3 and 5
Rscript scripts/ejs_simulations.R full results/generated/full

# Recreate Figures 1–3 and 5 from an existing full-run output directory
Rscript scripts/ejs_simulations.R figures results/generated/full

# Manuscript Figure 4: the ONLY official generator, using archived paired replications
Rscript scripts/figure4_article.R results/reference/replications.csv results/generated/figure4

# Table 3: independent degree-sensitivity experiment (m=3 vs m=4)
Rscript scripts/table3_degree_sensitivity.R results/generated/table3_degree_sensitivity

# Table 4: support-sensitivity independent audit (fast test or full)
Rscript scripts/table4_support_sensitivity.R test results/generated/table4_support_sensitivity/test
Rscript scripts/table4_support_sensitivity.R full results/generated/table4_support_sensitivity/full

# Local covariance alternatives: quick check or full experiment
Rscript scripts/ejs_simulations.R local_test results/generated/local_test
Rscript scripts/ejs_simulations.R local_full results/generated/local_full
```

The full and complementary Monte Carlo experiments can be time-consuming. The commands above **do not overwrite** the archived reference data in `results/reference/`.

### Outputs and manuscript mapping

| Manuscript item | Generator | Output |
| --- | --- | --- |
| Main simulation summaries and replication-level results (including Tables 1–2) | `ejs_simulations.R full` | `results/generated/full/summary.csv`, `replications.csv`, `adaptive_diagnostics.csv`, and covariance CSV files |
| Figures 1–3 | `ejs_simulations.R full` or `figures` | `results/generated/full/figure1.pdf`, `figure2.pdf`, `figure3.pdf` (also PNG files) |
| **Figure 4 (manuscript)** | **`figure4_article.R` only** | **`results/generated/figure4/figure4.pdf`**, `figure4.png`, `figure4_audit.csv` |
| Figure 5 | `ejs_simulations.R full` or `figures` | `results/generated/full/figure5.pdf` (also PNG) |
| Table 3 (degree sensitivity) | `table3_degree_sensitivity.R` | `results/generated/table3_degree_sensitivity/degree_power.csv`, `paired_power_differences.csv`, `null_size_check.csv`, `paired_rejections.csv` |
| Table 4 (support sensitivity) | `table4_support_sensitivity.R full` | `results/generated/table4_support_sensitivity/full/table4_audit.csv`, `table4_replications.csv` |
| Local alternatives | `ejs_simulations.R local_full` | `results/generated/local_full/local_power.csv`, `null_calibration.csv`, `local_power.pdf` |

**Figure 4 reference:** `results/reference/figure4.pdf` is the archived figure used for the manuscript. Regeneration writes to `results/generated/figure4/figure4.pdf`, so the reference is not overwritten. Its paired percentile bootstrap uses 2,000 resamples with seed `20261108`; compare the new figure visually and compare `figure4_audit.csv` with the archived replication data.

**Table 3:** the dedicated script reports the paired degree-3 versus degree-4 covariance-test sensitivity experiment; it does not produce a manuscript PDF table automatically. **Table 4:** the independent audit checks the published rounded MSEs statistically; it is not a bitwise replay of the historical random-number stream. For either table, compare the CSV values to the manuscript rather than assuming exact identity of independent Monte Carlo runs.

## Figure 4 conventions

- Panel A Y label: `Frequency of selecting the ratio estimator`.
- Panel B Y label: `MSE(theta_hat[n]^(ad,pr)) / MSE(theta_hat[n,3])` (rendered with R plotmath).
- Six sample sizes are equally spaced as categories: 250, 500, 1000, 2000, 4000, 10000.
- The 95% percentile intervals are computed using 2,000 resamples of Monte Carlo replications; the relative-MSE bootstrap is paired between estimators.
- The practical adaptive rule is a finite-sample implementation, not covered by the theoretical uniform guarantee.

## Important verification status

The repository contains the archived Monte Carlo results and reference figures associated with the manuscript. The author reports having executed and checked the R scripts prior to release. During preparation of this archive, the archived numerical outputs were also compared with the manuscript's four numerical tables; this comparison is not a substitute for an independent end-to-end R reproduction.

The archived manuscript Figure 4 is `results/reference/figure4.pdf`. To check the figure generator, run the dedicated command above and compare the newly generated figure with this reference. Bootstrap intervals may differ across R versions or random-number generator settings despite the fixed seed.

The practical adaptive selector and its calibration are empirical and do not inherit the theoretical adaptive estimator's uniform guarantees. The Table 4 support-sensitivity experiment is an independent Monte Carlo audit, not a bitwise recreation of the historical random-number stream. Independent end-to-end reproduction in a separate R environment has not been documented.

The manuscript itself is not included in this repository. The archived-results comparison does not constitute an exhaustive verification of all theoretical claims or proofs.
