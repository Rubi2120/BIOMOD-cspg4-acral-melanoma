# Computational reproducibility

## 1. Scope

This repository preserves the bioinformatic analysis developed for Pandora's Box, BIOMOD 2026.

The objective is to reproduce the statistical results and original English-language figures associated with the project, using the specified prepared transcriptomic and clinical inputs.

The workflow does not reproduce upstream sequencing, read alignment, expression quantification, or the original transcriptomic clustering procedure.

The computational verification was performed in November 2026.

## 2. Required environment

The validated workflow was executed using R 4.5.2.

The principal required R packages are:

- SummarizedExperiment
- ggplot2
- dplyr
- brunnermunzel
- survival

The package versions and platform information recorded during the integrated execution are provided in R_environment.txt.

The recorded environment describes the configuration used for verification. It is not an automatically provisioned or containerized software environment.

Changes in R versions, package versions, operating systems, graphics devices, or fonts may affect numerical or graphical reproducibility.

## 3. Required inputs

The workflow requires the 11 prepared input files described in DATA_AVAILABILITY.md.

The authorized input directory must remain separate from the public repository.

The required directory structure and original input filenames must be preserved.

The workflow does not automatically retrieve restricted datasets or reconstruct the prepared TCGA SummarizedExperiment object.

## 4. Executing the workflow

From the repository root, run:

```bash
Rscript --vanilla run_all.R \
  "/path/to/authorized/input_project" \
  "/path/to/private/new_analysis_run"
```

The second argument specifies a new output directory.

The output directory must not already exist, and its parent directory must exist.

The run directory must remain outside the public repository and the original input project.

The workflow creates its working structure and stages the required inputs and historical reference.

It then executes the 13 analysis scripts in the prescribed order.

The original input datasets are not intended to be modified.

## 5. Analysis sequence

| Script | Purpose |
| --- | --- |
| 00_prepare_tables.R | Prepare sample-level expression and tissue-reference tables |
| 00c_prepare_control_matrices.R | Harmonize expression matrices for internal controls |
| 00d_prepare_clinical_clusters.R | Reconstruct the 44-sample clinical and cluster annotations |
| 05_main_expression_analysis.R | Principal expression analysis, sensitivity analyses, and internal checks |
| 06_clinical_analysis.R | Clinical and molecular association analyses |
| 07_internal_expression_controls.R | Raw-expression and percentile-based internal controls |
| 08_normal_tissue_expression.R | Normal-tissue expression and vascular-marker correlations |
| 09_statistical_test_comparison.R | Comparison of alternative rank-based statistical tests |
| 10_TCGA_sensitivity_analysis.R | TCGA cohort-selection sensitivity analyses |
| 11_primary_vs_metastatic_acral.R | Primary versus lymph-node metastatic acral melanoma |
| 12_acral_survival_biomod.R | Historical BIOMOD survival analysis |
| 13_figure1_summary_statistics.R | Expression-threshold proportions and confidence intervals |
| 14_figures_english.R | Reproduction of the three original English-language figures |

The scripts are executed sequentially, with intermediate data generated before the downstream analyses that require them.

## 6. Random-number generation

The principal bootstrap procedures use 2,000 replicates and seed 20261022.

The verified statistical scripts retain explicit random-number generator settings, including the L'Ecuyer-CMRG generator with Inversion and Rejection settings.

The figure-generation script uses Mersenne-Twister, Inversion, and Rejection with seed 20261022.

The separate figure-generation setting preserves the random jitter coordinates of the original figures.

The statistical and graphical random-number generator settings must not be interchanged.

The script-level definitions are the authoritative implementation of these settings.

## 7. Verification of intermediate data

The English preparation scripts were verified against the original intermediate objects.

The reconstructed control expression matrices had the following dimensions:

| Group | Genes | Samples |
| --- | ---: | ---: |
| G1a | 56,156 | 604 |
| G1b | 56,156 | 701 |
| G2 | 60,616 | 103 |
| G3 | 60,666 | 44 |

The harmonized intersection contained 55,598 genes.

The reconstructed matrices matched the original values and dimension identifiers.

The clinical preparation script reproduced the 44-sample annotation table and the original cluster assignments.

The three transcriptomic clusters contained 14, 16, and 14 samples.

These checks were performed without publishing individual-level matrices or clinical annotations.

## 8. Integrated workflow verification

An initial integrated execution completed all 13 scripts.

Its outputs were compared with the established English-language reference results.

The comparison verified:

- All 20 scientific result tables.
- All three original PNG figures.
- All 16 internal checks evaluated during that execution.

One additional historical P-value comparison was absent from the newly generated audit table because its reference file had not been staged in the integrated working directory.

This discrepancy affected the number of audit checks, not the scientific results.

## 9. Historical P-value reconciliation

The original main expression analysis conditionally compared its calculated P values with those of a previous analysis.

The historical reference was recovered from the original project.

All six historical P values agreed with those produced by the integrated workflow.

The comparison used all.equal() on log10-transformed P values with a tolerance of 1e-6.

The aggregate reference is provided in reference_previous_p_values.csv.

Although the historical audit label refers to Brunner-Munzel P values, the original comparison includes all six principal contrasts: five independent-group Brunner-Munzel tests and one paired Wilcoxon signed-rank test.

The main expression script was subsequently executed with the historical reference staged in the expected working location.

This additional execution reproduced all eight aggregate tables generated by the main expression script, including the complete audit table.

All 17 internal checks passed.

The final run_all.R incorporates staging of the historical reference and requires the complete 17-check audit.

Its syntax was verified after this modification.

The modified runner was not independently rerun from beginning to end after the reconciliation.

Therefore, the final verification combines the original complete integrated execution with the separately verified historical-reference reconciliation.

## 10. Final verification results

| Verification item | Result |
| --- | --- |
| Analysis scripts completed in integrated execution | 13/13 |
| Scientific result tables reproduced | 20/20 |
| Reconciled internal audit table reproduced | 1/1 |
| Total aggregate tables reproduced | 21/21 |
| Original PNG figures reproduced | 3/3 |
| Internal checks passed after reconciliation | 17/17 |
| Historical P values verified | 6/6 |
| Combined output checks passed | 24/24 |

Aggregate CSV tables were compared after importing them into R using all.equal() with a numerical tolerance of 1e-10.

The three PNG figures were verified using MD5 checksums and reproduced byte for byte in the validated environment.

Numerical equivalence of the CSV tables does not necessarily imply byte-for-byte identity of their text serialization.

The final reconciled audit records the original integrated execution as the source of the scientific table and figure checks, and the subsequent main-script execution as the source of the complete audit-table verification.

## 11. Preservation of the original BIOMOD analysis

The repository preserves the analytical implementation underlying the original BIOMOD results.

The historical survival model is retained without substitution by a later delayed-entry analysis.

The expanded TCGA sensitivity analysis is distinguished from the earlier delivered sensitivity results in TCGA_sensitivity_versions.md.

The original figure files were not modified to accommodate the English-language script reorganization.

The objective of these decisions is reproducibility of the reported work, not retrospective replacement of the original analytical choices.

## 12. Private outputs and public release

The workflow generates sample-level tables, intermediate expression matrices, clinical annotations, and diagnostic records in the private run directory.

These materials are not intended for public redistribution.

The public repository contains analysis code, aggregate reference tables, and documentation, with figure redistribution subject to the appropriate authorization.

The first and third original figures contain individual plotted expression observations.

Reproducibility verification does not establish permission to redistribute patient-level data or individual plotted values.

See DATA_AVAILABILITY.md for the data-access requirements and privacy considerations.


## 13. Additional statistical dependency

The validated environment includes FSA version 0.10.1.
Script 05 uses FSA::dunnTest conditionally for Dunn's sensitivity analysis.
FSA should be available to reproduce the validated computational branch.
The associated dunn.test package was recorded at version 1.4.2.
