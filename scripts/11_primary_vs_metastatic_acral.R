# BIOMOD 2026 - Pandora's Box
# CSPG4 expression in primary acral melanoma and
# lymph node metastases.
#
# Compare 44 and 65 primary tumors with 15 lymph node
# metastases. Exclude lung and skin metastases and
# recurrent tumors from these comparisons.
#
# Some patients contribute primary and metastatic samples.
# The main comparisons are therefore exploratory and
# include potentially non-independent observations.
#
# An additional sensitivity analysis excludes lymph node
# metastases linked to patients with primary samples.
#
# Patient linkage is inferred from identifier prefixes.
# This is validated for samples in the supplementary table,
# but remains an assumption for other samples.
#
# Individual-level results must never be publicly exported.

local({

  input_root <- getOption(
    "biomod.metastases.input_root", "."
  )

  output_dir <- getOption(
    "biomod.metastases.output_dir", "tables"
  )

  private_dir <- getOption(
    "biomod.metastases.private_dir",
    file.path(tempdir(), "biomod_metastases_private")
  )

  if (!requireNamespace("brunnermunzel", quietly = TRUE)) {
    stop("Required package: brunnermunzel")
  }

  # Restore the original BIOMOD random number configuration.

  RNGkind(
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    private_dir,
    recursive = TRUE,
    showWarnings = FALSE,
    mode = "0700"
  )

  a <- read.csv(
    file.path(
      input_root,
      "results_v2/tables/nivelA_CSPG4_TPMcrudo_86muestras.csv"
    ),
    check.names = FALSE
  )

  # Identify the sample-type column using the original
  # detection procedure and source vocabulary.

  type_column <- NULL

  for (cn in names(a)) {

    v <- as.character(a[[cn]])

    if (
      any(grepl("primar", v, ignore.case = TRUE)) &&
      any(grepl(
        "metast|lymph|ganglio|node",
        v,
        ignore.case = TRUE
      ))
    ) {

      type_column <- cn
      break
    }
  }

  if (is.null(type_column)) {
    stop("Sample-type column could not be identified.")
  }

  v <- tolower(as.character(a[[type_column]]))

  primary <- grepl("primar", v)

  lymph_node <- grepl(
    "lymph|ganglio|node", v
  ) & !primary

  cat(
    "Primary tumors:", sum(primary),
    "| Lymph node metastases:", sum(lymph_node),
    "| Primary tumors in main analysis:",
    sum(primary & a$en_44),
    "\n"
  )

  stopifnot(
    sum(primary) == 65,
    sum(lymph_node) == 15,
    sum(primary & a$en_44) == 44
  )

  x <- a$CSPG4_log2

  # Cliff's delta calculated using pairwise differences.

  cliff_delta <- function(p, q) {
    mean(sign(outer(p, q, "-")))
  }

  # Original Brunner-Munzel implementation.
  # Cross-check against the brunnermunzel package.

  bm <- function(p, q) {

    n1 <- length(p)
    n2 <- length(q)

    r <- rank(c(p, q))
    r1 <- rank(p)
    r2 <- rank(q)

    m1 <- mean(r[1:n1])
    m2 <- mean(r[-(1:n1)])

    v1 <- sum(
      (r[1:n1] - r1 - m1 + (n1 + 1) / 2)^2
    ) / (n1 - 1)

    v2 <- sum(
      (r[-(1:n1)] - r2 - m2 + (n2 + 1) / 2)^2
    ) / (n2 - 1)

    W <- n1 * n2 * (m2 - m1) /
      ((n1 + n2) * sqrt(n1 * v1 + n2 * v2))

    df <- (n1 * v1 + n2 * v2)^2 /
      (
        (n1 * v1)^2 / (n1 - 1) +
        (n2 * v2)^2 / (n2 - 1)
      )

    2 * pt(-abs(W), df)
  }

  set.seed(20261022)

  # Preserve the original bootstrap calculations and
  # the sequence of comparisons.

  compare <- function(label, p, q) {

    stopifnot(abs(
      bm(p, q) -
        brunnermunzel::brunnermunzel.test(p, q)$p.value
    ) < 1e-6)

    bs <- replicate(2000, {

      pp <- sample(p, replace = TRUE)
      qq <- sample(q, replace = TRUE)

      c(
        cliff_delta(pp, qq),
        median(pp) - median(qq)
      )
    })

    data.frame(
      contrast = label,
      n_primary = length(p),
      n_metastases = length(q),
      median_primary = median(p),
      median_metastases = median(q),
      cliff_delta = cliff_delta(p, q),
      CI95_lower = quantile(bs[1, ], .025),
      CI95_upper = quantile(bs[1, ], .975),
      p_raw = bm(p, q),
      row.names = NULL
    )
  }

  # Main comparisons: BH adjustment across two tests.

  results <- rbind(
    compare(
      "44 primary tumors vs lymph node metastases",
      x[primary & a$en_44],
      x[lymph_node]
    ),
    compare(
      "65 primary tumors vs lymph node metastases",
      x[primary],
      x[lymph_node]
    )
  )

  results$p_BH <- p.adjust(results$p_raw, "BH")

  write.csv(
    results,
    file.path(
      output_dir,
      "acral_primary_vs_lymph_node_metastases.csv"
    ),
    row.names = FALSE
  )

  print(results, digits = 6, row.names = FALSE)

  # Patient linkage and potential sample dependence.
  # The supplementary table validates the prefix-based
  # convention for matched samples.

  s <- read.csv(
    file.path(
      input_root,
      "data/raw/Supplementary_Table_1.csv"
    ),
    check.names = FALSE
  )

  xp <- merge(
    data.frame(
      Sample = sub("^PR", "PD", a$sample_id)
    ),
    s[, c("Sample", "Patient")],
    by = "Sample"
  )

  cat(
    "\nSamples matched to published supplementary data:",
    nrow(xp),
    "\n"
  )

  stopifnot(all(
    sub("[a-z]+$", "", xp$Sample) == xp$Patient
  ))

  a$Patient <- sub(
    "[a-z]+$",
    "",
    sub("^PR", "PD", a$sample_id)
  )

  a$type <- ifelse(
    primary,
    "primary",
    ifelse(lymph_node, "lymph_node", "other")
  )

  tt <- table(
    a$Patient,
    factor(
      a$type,
      levels = c("primary", "lymph_node", "other")
    )
  )

  overlap <- rownames(tt)[
    tt[, "primary"] > 0 &
    tt[, "lymph_node"] > 0
  ]

  cat(
    "Patients with primary and lymph node samples:",
    length(overlap),
    "\n"
  )

  # Private individual-level output.
  # Do not print or copy this table to the public repository.

  paired_samples <- a[
    a$Patient %in% overlap & a$type != "other",
    c(
      "Patient",
      "sample_id",
      "type",
      "en_44",
      "CSPG4_log2"
    )
  ]

  write.csv(
    paired_samples,
    file.path(
      private_dir,
      "primary_and_metastatic_patient_samples.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nPrimary sample counts per inferred patient:\n"
  )

  print(table(table(a$Patient[primary])))

  # Sensitivity analysis:
  # Exclude lymph node samples linked to patients who
  # also contributed primary tumor samples.
  #
  # This reduces known potential overlap between groups,
  # but does not establish complete sample independence.
  # P values remain unadjusted, as in the original analysis.

  metastases_no_overlap <- lymph_node &
    !(a$Patient %in% overlap)

  cat(
    "Lymph node metastases without identified overlap:",
    sum(metastases_no_overlap),
    "\n"
  )

  sensitivity <- rbind(
    compare(
      "44 primary tumors vs metastases without identified overlap",
      x[primary & a$en_44],
      x[metastases_no_overlap]
    ),
    compare(
      "65 primary tumors vs metastases without identified overlap",
      x[primary],
      x[metastases_no_overlap]
    )
  )

  print(sensitivity, digits = 6, row.names = FALSE)

  write.csv(
    sensitivity,
    file.path(
      output_dir,
      "acral_metastasis_overlap_sensitivity.csv"
    ),
    row.names = FALSE
  )

  stopifnot(
    nrow(results) == 2,
    nrow(sensitivity) == 2,
    length(overlap) == 5,
    sum(metastases_no_overlap) == 9
  )

  # Interpretation:
  # Negative Cliff's delta indicates lower expression in
  # primary tumors than in lymph node metastases.
  #
  # The small sample sizes, potential patient overlap and
  # unadjusted sensitivity analyses limit interpretation.
  # These results do not establish differences in protein
  # expression or metastatic behavior.

})
