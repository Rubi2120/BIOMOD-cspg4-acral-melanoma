# Bioinformatic and statistical methods

## 1. Study objective

The bioinformatic analysis of Pandora's Box was designed to characterize the expression of chondroitin sulfate proteoglycan 4 (CSPG4) in Mexican primary acral melanoma and evaluate its transcript-level distribution relative to normal skin, primary melanoma from TCGA-SKCM, and other normal tissue categories.

The analysis provides a computational rationale for investigating CSPG4 as a molecular target for an aptamer-functionalized nanovehicle intended for the localized delivery of ribociclib.

This study is a retrospective analysis of existing transcriptomic and clinical datasets. It does not establish protein-level expression, cell-surface accessibility, aptamer binding, internalization, or therapeutic efficacy.

## 2. Data sources and sample selection

Four groups were defined:

| Group | Description | Sample count |
| --- | --- | ---: |
| G1a | GTEx v8 non-sun-exposed suprapubic skin | 604 |
| G1b | GTEx v8 sun-exposed lower-leg skin | 701 |
| G2 | TCGA-SKCM primary melanoma | 103 |
| G3 | Mexican primary acral melanoma | 44 |

The GTEx skin groups included 470 donors with paired samples from both anatomical sites.

The TCGA-SKCM reference group was restricted to primary tumors. Metastatic specimens were excluded from the principal comparison.

The TCGA group was not assumed to represent an exclusively non-acral melanoma population.

The Mexican dataset originated from the study by Basurto-Lozada et al., Ancestry and somatic profile indicate acral melanoma origin and prognosis, Nature (2026).

The dataset supplied for the BIOMOD reanalysis contained 86 specimens, including 65 primary tumors, 19 metastatic specimens, and two recurrent tumors.

The main acral melanoma analysis used the 44 primary tumors associated with the quality-controlled transcriptomic clustering analysis. An additional sensitivity analysis included all 65 available primary specimens.

The original study collected specimens at the Instituto Nacional de Cancerologia in Mexico between 2017 and 2019.

RNA was extracted from formalin-fixed, paraffin-embedded material at the Wellcome Sanger Institute. The original RNA-seq workflow used the Illumina HiSeq 4000 platform, STAR version 2.5.0 against GRCh38, and HTSeq version 0.7.2.

No additional recruitment, specimen collection, RNA extraction, or sequencing was performed for the BIOMOD reanalysis.

## 3. Expression harmonization

CSPG4 was identified using the Ensembl gene identifier ENSG00000173546.

Ensembl version suffixes were removed when necessary to harmonize gene identifiers across datasets.

GTEx expression values originated from the v8 RNA-SeQC TPM matrices. TCGA expression values were extracted from the unstranded TPM assay of the prepared SummarizedExperiment object.

For the Mexican cohort, CSPG4 TPM values were extracted from the original 86-specimen expression matrix and restricted to the selected 44 primary tumors.

Expression was transformed as:

log2(TPM + 1)

The original TPM scale was retained for comparisons involving tissue-level median expression values and within-sample expression percentiles.

No cross-cohort batch correction was applied.

Consequently, differences between cohorts may reflect both biological variation and technical differences in specimen processing, sequencing, expression quantification, and cohort composition.

## 4. Principal statistical analysis

CSPG4 expression distributions were summarized using medians and interquartile ranges.

Six pairwise contrasts were evaluated among the four groups.

The five independent-group comparisons used the Brunner-Munzel test. Cliff's delta was calculated to quantify the direction and magnitude of distributional differences.

Negative Cliff's delta values indicate lower expression in the first group relative to the second group. Positive values indicate the opposite direction.

Independent-group effect-size confidence intervals were estimated using 2,000 bootstrap replicates.

The paired comparison between sun-exposed and non-sun-exposed GTEx skin used the Wilcoxon signed-rank test, restricted to the 470 matched donors.

The paired difference was summarized using the Hodges-Lehmann estimate and its confidence interval.

P values from the six principal contrasts were adjusted together using the Benjamini-Hochberg procedure.

An adjusted P value below 0.05 was used as the statistical significance threshold.

The comparisons were exploratory and should not be interpreted as prospective confirmatory tests.

## 5. Internal expression controls

Six genes were evaluated to contextualize CSPG4 expression across datasets:

- CSPG4: target gene.
- MLANA and PMEL: melanocytic lineage-associated genes.
- ACTB, GAPDH, and TBP: conventional reference genes.

The gene-level TPM matrices were harmonized using shared Ensembl identifiers.

The common expression background comprised 55,598 genes.

For each control gene, expression was examined using both log2(TPM + 1) and within-sample expression percentiles relative to the common gene background.

Percentile calculations used positive TPM values, with zero-expression genes assigned a percentile of zero.

Cliff's delta and Brunner-Munzel tests were used to compare the acral melanoma group with the two normal skin groups and the TCGA primary melanoma group.

Multiple-testing correction was applied separately to the 18 contrasts for each expression metric.

These analyses provide a sensitivity assessment of expression patterns but do not remove all cross-cohort technical variation.

## 6. Normal tissue expression

CSPG4 expression was examined across 54 GTEx tissue and cell-line categories.

Each GTEx category was represented by its median CSPG4 TPM value.

These medians were compared with the median CSPG4 expression of the 44 primary acral melanomas.

The analysis identified the number of normal tissue categories with median expression above the acral melanoma median.

Additional analyses estimated the proportion of acral melanomas exceeding selected normal tissue expression thresholds.

Exact binomial confidence intervals were used for these proportions.

The pooled normal skin reference was also used to calculate the 95th percentile of CSPG4 expression.

To contextualize the normal-tissue expression pattern, Spearman correlations were calculated between tissue-level CSPG4 expression and the vascular or mural cell-associated genes PDGFRB, RGS5, ACTA2, and MYH11.

The two cell-line categories were excluded from these correlations, leaving 52 tissue categories.

These tissue-level associations do not establish the cellular origin of CSPG4 expression.

## 7. Clinical and molecular associations

Clinical associations were evaluated within the 44-sample primary acral melanoma subset.

Existing transcriptomic cluster assignments were obtained from the processed expression and classification files associated with the original cohort analysis.

According to correspondence with the original analyst, the expression matrix had been processed using the cola adjust_matrix() function to remove low-variance genes before clustering.

The exact settings applied to the supplied processed matrix were not independently confirmed.

The existing cluster assignments were retained; de novo clustering was not performed for the BIOMOD reanalysis.

The clinical analysis used CSPG4 expression from the original TPM matrix.

Clinical and molecular variables included transcriptomic cluster, pathological stage, ulceration, and sex.

The following tests were applied:

- Kruskal-Wallis test for differences among transcriptomic clusters.
- Spearman correlation for the association with ordinal tumor stage.
- Brunner-Munzel tests for binary clinical variables.
- Cliff's delta for binary-group effect sizes.

The four principal clinical association P values were adjusted together using the Benjamini-Hochberg procedure.

Pairwise cluster comparisons were corrected separately.

The public descriptive table suppresses selected summary statistics for categories with fewer than five observations.

This suppression is a precautionary disclosure-control measure rather than proof of redistribution authorization.

## 8. Tumor purity and deconvolution

Spearman correlations were used to investigate associations between CSPG4 expression and available tumor purity estimates.

Purity estimates were available for 22 primary acral melanoma samples.

Existing deconvolution estimates were examined for the 44 primary acral melanoma samples.

The deconvolution analysis evaluated eight estimated cellular fractions.

These analyses are exploratory and do not directly identify which cells express CSPG4 within the tumors.

## 9. Statistical and cohort sensitivity analyses

The principal findings were evaluated using complementary statistical and cohort-selection approaches.

A Kruskal-Wallis omnibus test followed by Dunn's pairwise comparisons was included as an alternative rank-based analysis.

Wilcoxon rank-sum tests and Brunner-Munzel tests were also compared for independent-group contrasts.

The primary acral melanoma analysis was repeated using the broader set of 65 primary specimens.

TCGA sensitivity analyses evaluated the complete 103-sample primary melanoma group, a 92-sample skin-primary subset, and a 102-sample subset excluding the specifically identified TCGA-ER-A196 case.

The distinction between the original delivered TCGA sensitivity results and the expanded sensitivity table is documented in TCGA_sensitivity_versions.md.

These analyses were used to assess the stability of the principal expression patterns rather than to establish population-level diagnostic or predictive performance.

## 10. Primary versus metastatic acral melanoma

CSPG4 expression was compared between primary acral melanoma specimens and 15 lymph-node metastatic specimens.

The analysis included comparisons using both the 44-sample primary subset and the broader 65-sample primary dataset.

Brunner-Munzel tests and Cliff's delta were used to characterize differences between primary and metastatic specimens.

Potential overlap between patients contributing primary and metastatic specimens was investigated using the available sample annotations and specimen identifier patterns.

A sensitivity analysis excluded lymph-node metastatic specimens with identified overlap with the primary tumor samples.

The available annotations did not permit definitive independent validation of all patient-level relationships.

Patient-level linkage tables are retained privately.

## 11. Exploratory survival analysis

The historical BIOMOD survival analysis evaluated the association between continuous CSPG4 expression and overall survival in the primary acral melanoma cohort.

The model used unadjusted Cox proportional hazards regression, with time measured from diagnosis.

The analyzed cohort included 44 samples and 10 observed events.

A median-based expression grouping was used for a complementary log-rank comparison.

The proportional hazards assumption was assessed using Schoenfeld residuals.

The analysis reproduces the historical BIOMOD implementation, which did not incorporate delayed entry.

The small number of events, lack of multivariable adjustment, and potential limitations associated with the time origin preclude interpretation as independent prognostic validation.

A subsequent delayed-entry sensitivity analysis was not substituted for the original BIOMOD result.

## 12. Figure generation

Three principal figures were reproduced:

1. CSPG4 expression across normal skin, TCGA primary melanoma, and Mexican primary acral melanoma.
2. CSPG4 median expression across 54 GTEx categories relative to the acral melanoma median.
3. CSPG4 expression across selected clinical and molecular characteristics of the 44 primary acral melanomas.

The original figure dimensions, plotting parameters, and graphical settings were retained.

The first and third figures contain individual expression observations and require appropriate redistribution authorization before public release.

## 13. Computational reproducibility

The complete workflow consists of 13 R scripts orchestrated by run_all.R.

The workflow begins with previously quantified transcriptomic data and prepared sample annotations.

The statistical analyses use fixed random seeds and explicit random-number generator settings to preserve the verified results.

The principal bootstrap analyses use 2,000 replicates and seed 20261022.

The original integrated workflow reproduced 20 scientific result tables and all three original figures.

A subsequent execution of the main expression analysis restored the historical P-value comparison and reproduced the complete 17-check audit table.

The combined verification therefore covers 21 aggregate tables and three figures.

Detailed software versions, random-number generator settings, and verification procedures are documented in REPRODUCIBILITY.md.

## 14. Interpretation and limitations

The observed differences in CSPG4 transcript abundance support further investigation of CSPG4 in acral melanoma.

However, expression in multiple normal tissues demonstrates that CSPG4 is not a tumor-exclusive transcript.

Bulk RNA-seq measurements cannot establish protein abundance, membrane localization, epitope accessibility, aptamer binding, or cellular internalization.

The datasets differ in anatomical sampling, population composition, tissue processing, and expression quantification.

The normal skin references are suprapubic and lower-leg skin rather than anatomically matched acral skin.

Clinical, metastatic, and survival analyses are exploratory and limited by sample size, event counts, and available metadata.

The proposed localized ribociclib delivery strategy requires independent experimental validation of target accessibility, aptamer specificity, nanovehicle uptake, drug activity, and safety.
