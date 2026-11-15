# Data availability and input requirements

## Overview

The Pandora's Box bioinformatic analysis integrates previously generated transcriptomic and clinical data from GTEx v8, TCGA-SKCM, and a Mexican acral melanoma cohort.

The repository provides analysis scripts, aggregate results, and reproducibility documentation. Individual-level transcriptomic matrices, clinical records, sample identifiers, and intermediate objects are not included.

## 1. GTEx v8

The Genotype-Tissue Expression (GTEx) Project provides the normal-tissue expression reference.

The primary comparison includes:

- 604 non-sun-exposed suprapubic skin samples.
- 701 sun-exposed lower-leg skin samples.
- 470 donors with paired samples from both skin categories.

The normal-tissue analysis additionally uses median expression values across 54 GTEx categories, including 52 tissue categories and two cell-line categories.

Expression values were obtained from the GTEx v8 RNA-SeQC gene-level TPM matrices, based on GRCh38 and GENCODE v26.

Resources:

- GTEx Portal: https://gtexportal.org/
- GTEx Consortium, Science (2020): https://doi.org/10.1126/science.aaz1776

The two GTEx expression files are distributed through the GTEx Portal. The skin sample metadata file used by this workflow is a prepared input containing the sample annotations required for the analysis.

## 2. TCGA-SKCM

The Cancer Genome Atlas Skin Cutaneous Melanoma project provides the reference group of 103 primary melanoma specimens.

Metastatic specimens were excluded from the primary expression comparison. The analysis uses gene-level unstranded TPM values from a previously prepared SummarizedExperiment object.

The TCGA reference group should not be interpreted as an exclusively non-acral melanoma cohort. Sensitivity analyses examine alternative sample selections.

Resources:

- GDC Data Portal: https://portal.gdc.cancer.gov/
- TCGA melanoma study: https://doi.org/10.1016/j.cell.2015.05.044
- GDC mRNA expression pipeline: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/

The prepared TCGA object is a required input. This repository does not reconstruct that object from sequencing reads or automatically download it from GDC.

## 3. Mexican acral melanoma cohort

The Mexican cohort originates from the study:

Basurto-Lozada et al. Ancestry and somatic profile indicate acral melanoma origin and prognosis. Nature (2026).

https://doi.org/10.1038/s41586-025-09967-z

The original study describes specimens collected at the Instituto Nacional de Cancerologia in Mexico between 2017 and 2019, with informed consent.

RNA was extracted from formalin-fixed, paraffin-embedded specimens. The original investigators performed RNA sequencing and expression quantification.

The expression dataset received for the BIOMOD reanalysis contained 86 specimens:

- 65 primary tumors.
- 19 metastatic specimens.
- 2 recurrent tumors.

The primary BIOMOD analysis uses a subset of 44 primary tumors associated with the quality-controlled transcriptomic clustering analysis. A sensitivity analysis includes all 65 available primary specimens.

The BIOMOD expression analysis uses the original TPM values for CSPG4 in the selected samples rather than the adjusted expression values used for transcriptomic clustering.

Processed expression matrices and existing cluster assignments are used to reconstruct the clinical and molecular annotations of the 44-sample subset.

BIOMOD did not recruit new patients, collect additional specimens, or perform new RNA sequencing for this reanalysis.

The cohort files supplied to the team are not redistributed through this repository. Access must be obtained through the appropriate data custodians and applicable authorization procedures.

## 4. Required input files

The complete workflow requires the following 11 files inside an authorized input project directory.

The original filenames are retained for compatibility with the verified analysis scripts.

| Relative path | Description |
| --- | --- |
| data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz | GTEx v8 sample-level gene TPM matrix |
| data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz | GTEx v8 median gene expression by tissue category |
| data/GTEx_skin_samples_v8.csv | Prepared GTEx skin sample metadata |
| data/raw/GDCdata/tcga_skcm_se_preparado.rds | Prepared TCGA-SKCM SummarizedExperiment object |
| data/raw/acral_MX_TPMs_qcpassed.tsv | Original Mexican cohort TPM matrix; 86 specimens |
| data/raw/acral_MX_metadata_.csv | Mexican cohort sample metadata |
| data/raw/acralMX_expression_data_update44s.rds | Processed expression matrix for the 44-sample subset |
| data/raw/acralMX_metadata_cc_update.csv | Clinical metadata for the 44-sample subset |
| data/raw/classification_cc_SDpam.rds | Existing transcriptomic cluster assignments |
| data/raw/Supplementary_Table_1.csv | Cohort supplementary clinical annotations |
| data/raw/Sample_pass_deconvolution_data.csv | Existing sample-level deconvolution estimates |

The file acral_MX_TPMs_qcpassed.tsv is read using its original comma-separated format despite the filename extension.

The exact prepared inputs are necessary for numerical reproduction of the reported analysis. Alternative downloads, updated annotations, or independently reconstructed matrices may produce different results.

## 5. Reproduction and output handling

The input project must remain separate from the public repository.

Run the workflow using:

```bash
Rscript --vanilla run_all.R \
  "/path/to/authorized/input_project" \
  "/path/to/private/new_analysis_run"
```

The output directory must not already exist.

The workflow prepares intermediate tables, expression matrices, and clinical annotations before executing the statistical analyses and figure-generation scripts.

Individual-level outputs and intermediate objects must remain in the private analysis environment.

The aggregate tables and figures are generated separately from these private intermediates.

## 6. Privacy and redistribution

The following materials are not intended for inclusion in the public GitHub repository:

- Original or processed patient-level transcriptomic matrices.
- Clinical metadata and sample identifiers.
- Individual survival records.
- Patient-level tumor purity and deconvolution estimates.
- Sample-level expression tables.
- Intermediate RDS and RData objects containing individual-level data.
- Local package libraries and private working directories.

The public availability of an associated scientific publication does not automatically establish permission to redistribute the particular data files supplied to the BIOMOD team.

Aggregate summaries must also be reviewed for disclosure risks, especially in clinical categories with very small sample counts.

The original expression comparison and clinical figures contain individual plotted observations. Their redistribution must be evaluated under the applicable data-sharing authorization before public release.

## 7. Scope of reproducibility

The repository documents the BIOMOD reanalysis of existing quantified expression data.

Reproducibility refers to the computational workflow beginning with the specified prepared inputs. It does not imply independent reproduction of the original cohort recruitment, laboratory procedures, sequencing, read alignment, or transcriptomic clustering.

The scripts and aggregate reference outputs are intended to support scientific review and computational verification, subject to authorized access to the necessary input datasets.
