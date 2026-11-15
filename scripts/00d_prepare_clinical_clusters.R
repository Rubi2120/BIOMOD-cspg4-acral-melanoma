# BIOMOD 2026 - Pandora's Box
# Reconstruct the original 44-sample clinical table.
#
# This script reproduces the original Level A preparation
# procedure using the supplied expression, classification
# and clinical metadata files.
#
# Cluster assignments are imported from the original
# classification object. Consensus clustering is not
# repeated or modified.
#
# IMPORTANT:
# The CSPG4 expression column in this intermediate table
# contains values from the processed 44-sample matrix.
#
# Script 06 uses original TPM values from the 86-sample
# matrix, restricted to the 44-sample subset, and imports
# only the cluster assignments from this table.
#
# The resulting table contains individual-level data
# and must remain private.

local({

  input_root <- getOption(
    "biomod.prepare.input_root", "."
  )

  work_root <- getOption(
    "biomod.prepare.work_root"
  )

  if (is.null(work_root)) {
    stop(
      "Set biomod.prepare.work_root ",
      "to a private directory."
    )
  }

  output_dir <- file.path(
    work_root,
    "results_v2/tables"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE,
    mode = "0700"
  )

  stopifnot(
    normalizePath(work_root) !=
      normalizePath(input_root)
  )

  output_file <- file.path(
    output_dir,
    "nivelA_CSPG4_44muestras.csv"
  )

  if (file.exists(output_file)) {
    stop(
      "Clinical intermediate table already exists. ",
      "Use a new private workspace."
    )
  }

  # -------------------------------------------------------
  # Load the original supplied files
  # -------------------------------------------------------

  expr44 <- readRDS(file.path(
    input_root,
    "data/raw/acralMX_expression_data_update44s.rds"
  ))

  meta44 <- read.csv(
    file.path(
      input_root,
      "data/raw/acralMX_metadata_cc_update.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  classification <- readRDS(file.path(
    input_root,
    "data/raw/classification_cc_SDpam.rds"
  ))

  gene_id <- "ENSG00000173546"

  # -------------------------------------------------------
  # Validate input structure and sample correspondence
  # -------------------------------------------------------

  stopifnot(
    is.matrix(expr44),
    is.numeric(expr44),
    gene_id %in% rownames(expr44),
    "sample_id" %in% names(meta44),
    is.data.frame(classification),
    "class" %in% names(classification),
    ncol(expr44) == 44,
    nrow(meta44) == 44,
    nrow(classification) == 44
  )

  stopifnot(
    !anyDuplicated(colnames(expr44)),
    !anyDuplicated(meta44$sample_id),
    !anyDuplicated(rownames(classification))
  )

  stopifnot(
    setequal(
      colnames(expr44),
      meta44$sample_id
    ),
    setequal(
      colnames(expr44),
      rownames(classification)
    )
  )

  stopifnot(
    !anyNA(expr44),
    all(is.finite(expr44)),
    all(expr44 >= 0)
  )

  # -------------------------------------------------------
  # Extract CSPG4 from the processed expression matrix.
  #
  # Preserve its original scale. Do not apply another
  # logarithmic transformation.
  # -------------------------------------------------------

  expression <- data.frame(
    sample_id = colnames(expr44),
    CSPG4_log2TPM1 = as.numeric(
      expr44[gene_id, ]
    ),
    stringsAsFactors = FALSE
  )

  # Import original cluster assignments using
  # classification row names as sample identifiers.

  clusters <- data.frame(
    sample_id = rownames(classification),
    cluster = classification$class,
    stringsAsFactors = FALSE
  )

  # Preserve the original merge procedure and ordering.

  results <- merge(
    expression,
    clusters,
    by = "sample_id",
    all = TRUE,
    sort = FALSE
  )

  results <- merge(
    results,
    meta44,
    by = "sample_id",
    all = TRUE,
    sort = FALSE
  )

  # -------------------------------------------------------
  # Verify the reconstructed table
  # -------------------------------------------------------

  expected_columns <- c(
    "sample_id",
    "CSPG4_log2TPM1",
    "cluster",
    "mutation_status",
    "stage",
    "sex",
    "age",
    "ulceration",
    "WGD",
    "TMB",
    "tumor_purity"
  )

  stopifnot(
    nrow(results) == 44,
    identical(
      names(results),
      expected_columns
    ),
    !anyDuplicated(results$sample_id),
    !anyNA(results$CSPG4_log2TPM1),
    !anyNA(results$cluster)
  )

  cluster_counts <- table(
    results$cluster
  )

  stopifnot(
    identical(
      names(cluster_counts),
      c("1", "2", "3")
    ),
    identical(
      as.integer(cluster_counts),
      c(14L, 16L, 14L)
    )
  )

  # Verify sample-by-sample assignment correspondence.

  expected_clusters <- classification$class[
    match(
      results$sample_id,
      rownames(classification)
    )
  ]

  stopifnot(identical(
    as.numeric(results$cluster),
    as.numeric(expected_clusters)
  ))

  # Write individual-level data only to the private
  # intermediate directory.

  write.csv(
    results,
    output_file,
    row.names = FALSE
  )

  Sys.chmod(
    output_file,
    mode = "0600"
  )

  cat(
    "\nOriginal cluster assignments imported successfully.\n"
  )

  cat(
    "Clinical table dimensions:",
    paste(dim(results), collapse = " x "),
    "\n"
  )

  cat("\nAggregate cluster counts:\n")

  print(cluster_counts)

  cat(
    "\nPrivate clinical intermediate saved to:\n",
    output_file,
    "\n"
  )

})
