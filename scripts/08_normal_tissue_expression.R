# BIOMOD 2026 - Pandora's Box
# CSPG4 expression across normal tissue categories.
#
# Objectives:
# 1. Estimate the proportion of primary acral melanomas
#    with CSPG4 expression above each GTEx tissue median.
# 2. Compare expression with the 95th percentile of
#    pooled non-acral skin samples.
# 3. Explore correlations between CSPG4 and vascular
#    marker expression across GTEx tissues.
#
# All analyses are descriptive and based on bulk RNA
# expression. They do not establish protein abundance,
# cellular specificity, target accessibility or safety.
#
# The individual-level expression matrices are restricted.
# Only aggregate results are exported.

local({

  suppressPackageStartupMessages(library(dplyr))

  input_root <- getOption(
    "biomod.tissues.input_root", "."
  )

  output_dir <- getOption(
    "biomod.tissues.output_dir", "tables"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # 1. Load expression data -----------------------------------

  cs <- "ENSG00000173546"

  input_file <- file.path(
    input_root,
    "results_v2/intermedios/ctl_tpm.rds"
  )

  if (!file.exists(input_file)) {
    stop("Missing expression matrices: ", input_file)
  }

  rr <- readRDS(input_file)

  tpm <- rr$tpm

  # Original TPM values for 44 primary acral melanomas.

  acral <- as.numeric(tpm$G3[cs, ])

  # Pool the two non-acral skin groups.

  skin <- c(
    as.numeric(tpm$G1a[cs, ]),
    as.numeric(tpm$G1b[cs, ])
  )

  # GTEx median expression matrix.

  md <- rr$med

  gene_column <- if ("Name" %in% names(md)) {
    "Name"
  } else {
    names(md)[1]
  }

  gene_row <- which(
    sub("\\..*$", "", md[[gene_column]]) == cs
  )

  stopifnot(length(gene_row) == 1)

  tissue_columns <- setdiff(
    names(md),
    c(gene_column, "Description")
  )

  tissue_medians <- unlist(
    md[gene_row, tissue_columns]
  )

  stopifnot(
    length(tissue_medians) == 54,
    length(acral) == 44,
    length(skin) == 1305
  )

  # Consistency check against the original analysis.
  # The highest GTEx median is approximately 246.2 TPM
  # in tibial artery.

  if (abs(max(tissue_medians) - 246.2) > 1) {

    warning(
      "Maximum tissue median is ",
      round(max(tissue_medians), 1),
      " TPM; expected approximately 246.2 TPM."
    )
  }

  # 2. Normal tissue expression comparisons ------------------
  #
  # For each GTEx category:
  # - Median CSPG4 expression in TPM.
  # - Number and percentage of acral tumors above that median.
  # - Exact binomial 95% confidence interval.
  # - Ratio of the acral median to the GTEx category median.

  tissue_summary <- bind_rows(lapply(
    names(tissue_medians),
    function(t) {

      u <- tissue_medians[[t]]

      k <- sum(acral > u)

      n <- length(acral)

      ci <- binom.test(k, n)$conf.int

      data.frame(
        tissue = t,
        cell_line = grepl("^Cells", t),
        median_TPM = u,
        k = k,
        n = n,
        percent_acral_above_median = 100 * k / n,
        CI95_lower_percent = 100 * ci[1],
        CI95_upper_percent = 100 * ci[2],
        acral_to_tissue_median_ratio = median(acral) / u
      )
    }
  )) |>
    arrange(desc(median_TPM))

  stopifnot(
    nrow(tissue_summary) == 54,
    sum(tissue_summary$cell_line) == 2,
    sum(tissue_summary$median_TPM > median(acral)) == 10
  )

  print(as.data.frame(head(tissue_summary, 15)))

  write.csv(
    tissue_summary,
    file.path(
      output_dir,
      "normal_tissue_expression_summary.csv"
    ),
    row.names = FALSE
  )

  # Reference comparison with pooled non-acral skin.
  # This is a descriptive expression threshold, not a
  # validated threshold for therapeutic targeting.

  skin_p95 <- unname(quantile(skin, .95))

  k95 <- sum(acral > skin_p95)

  message(sprintf(
    paste0(
      "Reference: %d of %d acral tumors (%.1f%%) ",
      "exceed the pooled skin 95th percentile (%.2f TPM)."
    ),
    k95,
    length(acral),
    100 * k95 / length(acral),
    skin_p95
  ))

  message(
    "GTEx cell-line categories: ",
    paste(
      tissue_summary$tissue[tissue_summary$cell_line],
      collapse = ", "
    )
  )

  # 3. Exploratory vascular marker correlations --------------
  #
  # Evaluate whether tissue-level CSPG4 expression
  # covaries with pericyte and smooth muscle markers.

  markers <- c(
    CSPG4 = cs,
    PDGFRB = "ENSG00000113721",
    RGS5 = "ENSG00000143248",
    ACTA2 = "ENSG00000107796",
    MYH11 = "ENSG00000133392"
  )

  positions <- match(
    markers,
    sub("\\..*$", "", md[[gene_column]])
  )

  if (any(is.na(positions))) {
    stop("One or more marker genes are missing.")
  }

  # Verify that gene symbols correspond to Ensembl IDs.

  if (
    "Description" %in% names(md) &&
    !identical(
      unname(as.character(md$Description[positions])),
      names(markers)
    )
  ) {

    stop(
      "Marker gene symbols do not match the Ensembl IDs: ",
      paste(
        as.character(md$Description[positions]),
        collapse = ", "
      )
    )
  }

  # Exclude the two cell-line categories.
  # Retain 52 GTEx tissues.

  tissue_only <- tissue_columns[
    !grepl("^Cells", tissue_columns)
  ]

  M <- t(as.matrix(
    md[positions, tissue_only]
  ))

  colnames(M) <- names(markers)

  vascular_markers <- setdiff(
    names(markers),
    "CSPG4"
  )

  vascular_correlations <- data.frame(
    marker = vascular_markers,
    spearman_rho = vapply(
      vascular_markers,
      function(k) {
        cor(
          M[, "CSPG4"],
          M[, k],
          method = "spearman"
        )
      },
      numeric(1)
    ),
    n_tissues = nrow(M)
  )

  stopifnot(
    nrow(M) == 52,
    nrow(vascular_correlations) == 4
  )

  print(vascular_correlations)

  write.csv(
    vascular_correlations,
    file.path(
      output_dir,
      "vascular_marker_correlations.csv"
    ),
    row.names = FALSE
  )

  # Interpretation:
  #
  # Tissue-level correlations are compatible with a
  # vascular contribution to CSPG4 expression, but they
  # cannot identify the cells expressing the transcript.
  #
  # No P values or confidence intervals are calculated
  # for these correlations. Tissue medians are not
  # independent biological observations.
  #
  # Normal tissue expression does not establish the
  # accessibility of CSPG4 to an aptamer-functionalized
  # nanovehicle or predict off-target toxicity.

  cat(
    "\nGTEx categories:",
    nrow(tissue_summary),
    "\nTissues used for correlations:",
    nrow(M),
    "\nCategories above the acral median:",
    sum(tissue_summary$median_TPM > median(acral)),
    "\n"
  )

  message(
    "Aggregate tissue results saved to: ",
    output_dir
  )

})
