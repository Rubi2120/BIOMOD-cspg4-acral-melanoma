# Pandora’s Box BIOMOD 2026

## CSPG4 expression analysis in Mexican acral melanoma

This repository contains the bioinformatic analysis developed for **Pandora’s Box**, a BIOMOD 2026 project investigating an aptamer-functionalized nanovehicle targeting chondroitin sulfate proteoglycan 4 (CSPG4) for the proposed localized delivery of ribociclib in acral melanoma.

The analysis characterizes CSPG4 transcript abundance in Mexican primary acral melanoma relative to non-acral normal skin and primary melanoma from TCGA-SKCM. Complementary analyses examine expression across normal tissue categories, internal expression controls, clinical characteristics, and exploratory sensitivity analyses.

**Scope:** This is a retrospective reanalysis of existing transcriptomic and clinical data. It does not establish protein abundance, cell-surface accessibility, aptamer binding or internalization, therapeutic selectivity, or treatment efficacy.

## Data and study design

| Group | Dataset | Samples |
| --- | --- | ---: |
| G1a | GTEx v8 non-sun-exposed suprapubic skin | 604 |
| G1b | GTEx v8 sun-exposed lower-leg skin | 701 |
| G2 | TCGA-SKCM primary melanoma | 103 |
| G3 | Mexican primary acral melanoma | 44 |

The two GTEx skin groups include 470 paired donors. The main acral analysis uses 44 primary tumors from the supplied 86-sample Mexican cohort; a sensitivity analysis includes all 65 primary tumor samples. G2 is the primary-tumor subset of TCGA-SKCM and is not restricted exclusively to non-acral melanoma.

CSPG4 expression (Ensembl `ENSG00000173546`) is analyzed on the log2(TPM + 1) scale, using the original TPM values. No additional cross-study batch-effect correction is applied. The main comparisons use paired Wilcoxon signed-rank tests for paired skin samples and Brunner–Munzel tests for independent groups, with Hodges–Lehmann estimates or Cliff’s delta as appropriate. Multiple testing is addressed using the Benjamini–Hochberg procedure. The analyses are exploratory.

Study source and disease condition are confounded in the principal cross-cohort comparisons. Observed differences cannot therefore be attributed exclusively to biological variation.

## Main findings

| Group | Median CSPG4 expression, log2(TPM + 1) |
| --- | ---: |
| G1a | 3.576 |
| G1b | 3.830 |
| G2 | 7.673 |
| G3 | 5.511 |

CSPG4 transcript abundance was higher in the acral cohort than in either normal-skin group, but lower than in the TCGA primary-melanoma group. The corresponding Cliff’s delta estimates (first group relative to second) were −0.660 for G1a vs G3, −0.605 for G1b vs G3, and +0.479 for G2 vs G3.

The TCGA-versus-acral comparison is sensitive to the statistical test used: the adjusted Brunner-Munzel P value was approximately 1.5e-7, whereas the adjusted Dunn P value was 0.077. This difference is documented in the statistical sensitivity analyses.

In the acral cohort, 36 of 44 tumors (81.8%) exceeded the non-sun-exposed skin 75th percentile, while 34 of 44 (77.3%) exceeded the sun-exposed skin 75th percentile. Ten of the 54 GTEx tissue and cell-line categories had median CSPG4 expression above the acral tumor median. These findings indicate heterogeneous expression and expression in multiple normal tissues, supporting further investigation rather than demonstrating tumor-specific targeting.

## Repository structure

```text
.
├── README.md
├── run_all.R
├── scripts/       # Data preparation, statistical analyses and figures
├── tables/        # Aggregate results and analysis checks
├── figures/       # Authorized English-language project figures
└── docs/          # Methodological and reproducibility documentation
```

### Analysis workflow

| Script | Purpose |
| --- | --- |
| `00_prepare_tables.R` | Reconstruct four expression tables from the supplied datasets. |
| `00c_prepare_control_matrices.R` | Reconstruct complete TPM matrices for the four groups. |
| `00d_prepare_clinical_clusters.R` | Reconstruct the clinical intermediate table using supplied cluster assignments. |
| `05_main_expression_analysis.R` | Main CSPG4 comparisons, sensitivity and complementary analyses. |
| `06_clinical_analysis.R` | Clinical and molecular associations in acral melanoma. |
| `07_internal_expression_controls.R` | Melanocytic and reference-gene expression controls. |
| `08_normal_tissue_expression.R` | GTEx tissue medians and vascular-marker correlations. |
| `09_statistical_test_comparison.R` | Comparison of Dunn, Wilcoxon and Brunner–Munzel tests. |
| `10_TCGA_sensitivity_analysis.R` | Sensitivity to the composition of the TCGA cohort. |
| `11_primary_vs_metastatic_acral.R` | Primary tumors versus lymph-node metastases. |
| `12_acral_survival_biomod.R` | Original exploratory BIOMOD survival analysis. |
| `13_figure1_summary_statistics.R` | Counts, percentages and binomial confidence intervals. |
| `14_figures_english.R` | Regenerate the three English figures. |

## Reproducing the analysis

### Requirements

The workflow was verified using R 4.5.2 on macOS. Required R packages are `SummarizedExperiment`, `ggplot2`, `dplyr`, `brunnermunzel` and `survival`.

FSA version 0.10.1 is additionally required to reproduce the validated Dunn sensitivity-analysis branch. See docs/R_environment.txt for the complete recorded environment.

Install the required packages in a compatible R environment if they are not already available:

```r
install.packages(c("ggplot2", "dplyr", "brunnermunzel", "survival", "FSA"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SummarizedExperiment")
```

For exact reproduction, use the validated package versions and graphics environment. The workflow does not install packages automatically. It can also use an existing `local_R_library/` directory in the external input project.

### Run the complete workflow

Clone the repository and run the following command from its root, replacing both example paths:

```bash
Rscript --vanilla run_all.R \
  "/path/to/authorized/input_project" \
  "/path/to/private/new_analysis_run"
```

The input project must contain the required `data/` directory and original quantified expression matrices, sample metadata, classification object and supplementary files. The input project and private run directory must be outside the repository. The parent of the new run directory must exist, and the run directory itself must not already exist.

The program executes all 13 scripts in order. It creates a private workspace, links the supplied data directory without copying large matrices, and writes intermediate files and final outputs to the new run directory. The symbolic link is not a read-only filesystem restriction; the supplied analysis scripts are designed to read, not modify, the original data.

### Outputs and validation

The verified workflow produces 21 CSV tables (20 scientific result tables and one audit table) and three PNG figures. Individual-level intermediates and session information remain in the private run directory.

The integrated execution reproduced all 20 scientific tables and all three figures. A separate rerun of script 05 with the frozen historical P-value reference reproduced the 17-row audit table and its other seven outputs. Together, these checks verified all 24 reference outputs. The final runner was updated to stage the historical reference and require all 17 checks; its syntax was validated after that change, without repeating the entire integrated run.

Bootstrap analyses use the verified L’Ecuyer-CMRG random-number generator and seed `20261022`; figure generation uses Mersenne-Twister with the same seed. All three English PNGs were reproduced byte for byte using R 4.5.2 and ggplot2 4.0.3. Different graphics devices, fonts or software versions may affect exact image rendering.

The frozen six-value historical regression reference is stored in [`docs/reference_previous_p_values.csv`](docs/reference_previous_p_values.csv). The TCGA sensitivity analysis versions are described in [`docs/TCGA_sensitivity_versions.md`](docs/TCGA_sensitivity_versions.md).

## Data access and limitations

GTEx and TCGA inputs originate from their respective data resources. The Mexican acral melanoma cohort originates from the study by Basurto-Lozada and colleagues (Nature, 2026). The supplied processed 44-sample matrix and classification object are used to recover the original subset and cluster assignments; consensus clustering is not repeated.

The original cohort analysis code is available at:
https://github.com/CGBio-Lab/Mex-acral-exomes-transcriptomes

Access to the original study data is subject to the applicable repository terms, access controls and data-sharing authorizations. This repository is intended to distribute analysis code and authorized aggregate outputs, not individual-level clinical data, sample identifiers or complete expression matrices. The figures containing individual expression points, especially the clinical figure with small subgroups, require appropriate authorization for redistribution.

Cross-study differences in sample collection, tissue composition, library preparation, quantification and other technical factors limit biological interpretation. The normal skin samples are not anatomically matched to acral skin. Bulk RNA expression cannot identify the cellular source of CSPG4 or establish cell-surface protein availability. The exploratory survival analysis is unadjusted for clinical confounders and is not a validated prognostic model.

## Primary data resources

- [GTEx Portal](https://gtexportal.org/): GTEx v8 expression and tissue resources.
- [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/): TCGA-SKCM data.
- [Basurto-Lozada et al., Nature (2026)](https://doi.org/10.1038/s41586-025-09967-z): original acral melanoma study.
- [NCBI Gene: CSPG4](https://www.ncbi.nlm.nih.gov/gene/1464): gene annotation.

---

**Pandora’s Box, BIOMOD 2026**

This repository documents the bioinformatic component of the project and preserves the analytical results used for the BIOMOD presentation.

## Figure availability



The three original figures were reproduced during computational verification. The aggregate normal-tissue figure is included in this release. Figures containing individual-level observations are retained privately pending redistribution authorization.
