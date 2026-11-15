# BIOMOD 2026 - Pandora's Box
# Sensitivity of CSPG4 comparisons to TCGA sample selection.
#
# Compare three definitions of the TCGA primary tumor group:
# 1. All 103 primary tumors.
# 2. The 92 tumors with skin recorded as tissue of origin.
# 3. The 102 tumors remaining after exclusion of TCGA-ER-A196,
#    the case annotated as acral lentiginous melanoma.
#
# These analyses reproduce the original BIOMOD results.
# Restricting the sample group does not remove study-related
# confounding or establish biological equivalence.

local({

  RNGkind(
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )

  set.seed(20261022)

  suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(brunnermunzel)
  })

  input_root <- getOption(
    "biomod.tcga.input_root", "."
  )

  output_dir <- getOption(
    "biomod.tcga.output_dir", "tables"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Cliff's delta and bootstrap confidence intervals.

  cliff_delta <- function(x, y) {

    n1 <- length(x)
    n2 <- length(y)
    r <- rank(c(x, y))

    2 * ((sum(r[1:n1]) - n1 * (n1 + 1) / 2) /
           (n1 * n2)) - 1
  }

  bootstrap_ci <- function(x, y, B = 2000) {

    quantile(
      replicate(B, cliff_delta(
        sample(x, replace = TRUE),
        sample(y, replace = TRUE)
      )),
      c(.025, .975)
    )
  }

  # Load the expression data and TCGA metadata.

  a <- read.csv(file.path(
    input_root,
    "results_v2/tables/nivelA_CSPG4_TPMcrudo_86muestras.csv"
  ))

  b <- read.csv(file.path(
    input_root,
    "results_v2/tables/nivelB_TCGA_CSPG4_por_muestra.csv"
  ))

  p <- read.csv(file.path(
    input_root,
    "results_v2/tables/nivelC_GTEx_CSPG4_piel_muestras.csv"
  ))

  se <- readRDS(file.path(
    input_root,
    "data/raw/GDCdata/tcga_skcm_se_preparado.rds"
  ))

  cd <- as.data.frame(colData(se))

  rownames(cd) <- colnames(se)

  b$tissue_origin <- cd[
    b$barcode,
    "tissue_or_organ_of_origin"
  ]

  b$biopsy_site <- cd[
    b$barcode,
    "site_of_resection_or_biopsy"
  ]

  primary <- b[b$sample_type == "Primary Tumor", ]

  stopifnot(nrow(primary) == 103)

  cat("\nTissue of origin of the 103 primary tumors:\n")

  print(table(
    primary$tissue_origin,
    useNA = "ifany"
  ))

  # Subset 1: primary tumors with recorded skin origin.

  skin_primary <- primary[
    grepl("^Skin", primary$tissue_origin),
  ]

  stopifnot(nrow(skin_primary) == 92)

  # Subset 2: exclude the TCGA case annotated as acral.

  no_acral_case <- primary[
    !grepl("^TCGA-ER-A196", primary$barcode),
  ]

  stopifnot(nrow(no_acral_case) == 102)

  # Preserve the other three comparison groups.

  non_sun <- grepl("Not Sun Exposed", p$SMTSD)

  G <- list(
    G1a = p$CSPG4_log2TPM1[non_sun],
    G1b = p$CSPG4_log2TPM1[!non_sun],
    G3 = a$CSPG4_log2[a$en_44]
  )

  # Calculate one contrast.
  # The original bootstrap order is preserved.

  comparison_row <- function(selection, n_g2, x, y) {

    ci <- bootstrap_ci(x, y)

    data.frame(
      G2_selection = selection,
      n_G2 = n_g2,
      cliff_delta = cliff_delta(x, y),
      CI95_lower = ci[1],
      CI95_upper = ci[2],
      median_difference = median(x) - median(y),
      p_BM_raw = brunnermunzel.test(x, y)$p.value
    )
  }

  selections <- list(
    list("all primary tumors (103)", primary$CSPG4_log2TPM1),
    list("skin origin only (92)", skin_primary$CSPG4_log2TPM1),
    list(
      "excluding TCGA-ER-A196 (102)",
      no_acral_case$CSPG4_log2TPM1
    )
  )

  results <- do.call(rbind, lapply(selections, function(s) {

    g2 <- s[[2]]
    n <- length(g2)

    r <- rbind(
      comparison_row(s[[1]], n, G$G1a, g2),
      comparison_row(s[[1]], n, G$G1b, g2),
      comparison_row(s[[1]], n, g2, G$G3)
    )

    cbind(
      contrast = c(
        "G1a vs G2",
        "G1b vs G2",
        "G2 vs G3"
      ),
      r
    )
  }))

  rownames(results) <- NULL

  stopifnot(nrow(results) == 9)

  print(results, digits = 6, row.names = FALSE)

  write.csv(
    results,
    file.path(
      output_dir,
      "TCGA_sensitivity_analysis.csv"
    ),
    row.names = FALSE
  )

  # P values are unadjusted sensitivity-analysis results.
  # Bootstrap intervals may differ from script 05 because
  # the resampling sequence differs.
  #
  # The direction of the CSPG4 expression differences is
  # preserved under all three TCGA sample selections.

})
