# TCGA sensitivity analysis: version history

BIOMOD 2026 - Pandora's Box

## BIOMOD delivery version

The original delivery table
`sensibilidad_TCGA_solo_piel.csv`
contains six contrasts covering two TCGA selections:

- All 103 primary tumors.
- The 92 primary tumors with skin recorded as
  tissue or organ of origin.

Each selection includes three contrasts:
G1a vs G2, G1b vs G2, and G2 vs G3.

## Final analysis version

The final results table contains nine contrasts.
It retains the original six contrasts and adds
three contrasts using the 102 primary tumors
remaining after exclusion of TCGA-ER-A196,
the case annotated as acral lentiginous melanoma.

The English script
`scripts/10_TCGA_sensitivity_analysis.R`
reproduces this nine-row final version.

Its output is
`tables/TCGA_sensitivity_analysis.csv`.

## Numerical verification

The six delivery rows were matched to the final
results using contrast labels and TCGA group sizes.

All six rows matched exactly, including:

- Sample sizes.
- Cliff's delta estimates.
- Bootstrap 95% confidence intervals.
- Median differences.
- Unadjusted Brunner-Munzel P values.

The English nine-row table also reproduced the
complete original final results with a numerical
comparison tolerance of 1e-10.

The three additional contrasts are retained as
a complementary sensitivity analysis and do not
replace or modify the original delivery results.

## Interpretation

The direction of the CSPG4 expression differences
is preserved under all three TCGA selections.

These exploratory analyses do not eliminate
study-related confounding or establish differences
in CSPG4 protein expression or targeting specificity.
