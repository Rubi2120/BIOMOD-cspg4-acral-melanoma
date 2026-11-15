# BIOMOD 2026 - Pandora's Box
#
# Complete reproducible workflow for the CSPG4 analysis.
#
# Usage:
#
# Rscript --vanilla run_all.R INPUT_ROOT RUN_DIR
#
# INPUT_ROOT:
#   Directory containing the original data/ structure.
#
# RUN_DIR:
#   New private directory for all generated outputs.
#   This directory must not already exist.
#
# Optional third argument:
#   Repository root, if running a copy of this script
#   from a different location.
#
# The workflow starts from previously quantified TPM
# matrices and supplied clinical/classification data.
#
# It does not download restricted data, perform RNA
# sequencing, align reads or repeat consensus clustering.
#
# No individual-level data are distributed by this script.
# All generated outputs remain in the private run directory.

main <- function() {

  args <- commandArgs(trailingOnly = TRUE)

  if (!length(args) %in% c(2L, 3L)) {

    stop(
      "Usage: Rscript --vanilla run_all.R ",
      "INPUT_ROOT RUN_DIR [REPOSITORY_ROOT]"
    )
  }

  # Restrict permissions on newly created files.

  Sys.umask("0077")

  input_root <- normalizePath(
    path.expand(args[1]),
    mustWork = TRUE
  )

  requested_run <- path.expand(args[2])

  run_parent <- normalizePath(
    dirname(requested_run),
    mustWork = TRUE
  )

  run_root <- file.path(
    run_parent,
    basename(requested_run)
  )

  # Locate the repository independently of the current
  # working directory.

  if (length(args) == 3L) {

    repo <- normalizePath(
      path.expand(args[3]),
      mustWork = TRUE
    )

  } else {

    script_arg <- grep(
      "^--file=",
      commandArgs(trailingOnly = FALSE),
      value = TRUE
    )

    stopifnot(length(script_arg) == 1)

    script_file <- sub(
      "^--file=", "", script_arg
    )

    repo <- dirname(normalizePath(
      script_file,
      mustWork = TRUE
    ))
  }

  # Never create the private run directory inside the
  # public repository or original input project.

  inside_directory <- function(path, root) {

    identical(path, root) ||
      startsWith(path, paste0(root, "/"))
  }

  stopifnot(
    dir.exists(input_root),
    dir.exists(repo),
    !inside_directory(run_root, repo),
    !inside_directory(run_root, input_root),
    !inside_directory(input_root, repo)
  )

  if (file.exists(run_root)) {
    stop("RUN_DIR already exists. Choose a new directory.")
  }

  # -----------------------------------------------------
  # 1. Verify required input files.
  # -----------------------------------------------------

  required_inputs <- c(

    "data/raw/acral_MX_TPMs_qcpassed.tsv",

    "data/raw/acral_MX_metadata_.csv",

    "data/raw/acralMX_expression_data_update44s.rds",

    "data/raw/acralMX_metadata_cc_update.csv",

    "data/raw/classification_cc_SDpam.rds",

    "data/raw/GDCdata/tcga_skcm_se_preparado.rds",

    paste0(
      "data/raw/",
      "GTEx_Analysis_2017-06-05_v8_",
      "RNASeQCv1.1.9_gene_tpm.gct.gz"
    ),

    paste0(
      "data/raw/",
      "GTEx_Analysis_2017-06-05_v8_",
      "RNASeQCv1.1.9_gene_median_tpm.gct.gz"
    ),

    "data/GTEx_skin_samples_v8.csv",

    "data/raw/Supplementary_Table_1.csv",

    "data/raw/Sample_pass_deconvolution_data.csv"
  )

  missing_inputs <- required_inputs[
    !file.exists(file.path(
      input_root,
      required_inputs
    ))
  ]

  if (length(missing_inputs)) {

    stop(
      "Missing required input files:\n",
      paste(missing_inputs, collapse = "\n")
    )
  }

  # -----------------------------------------------------
  # 2. Verify all analysis scripts.
  # -----------------------------------------------------

  scripts <- c(
    "00_prepare_tables.R",
    "00c_prepare_control_matrices.R",
    "00d_prepare_clinical_clusters.R",
    "05_main_expression_analysis.R",
    "06_clinical_analysis.R",
    "07_internal_expression_controls.R",
    "08_normal_tissue_expression.R",
    "09_statistical_test_comparison.R",
    "10_TCGA_sensitivity_analysis.R",
    "11_primary_vs_metastatic_acral.R",
    "12_acral_survival_biomod.R",
    "13_figure1_summary_statistics.R",
    "14_figures_english.R"
  )

  script_paths <- file.path(
    repo, "scripts", scripts
  )

  missing_scripts <- scripts[
    !file.exists(script_paths)
  ]

  if (length(missing_scripts)) {

    stop(
      "Missing scripts:\n",
      paste(missing_scripts, collapse = "\n")
    )
  }

  # Use an existing local R library when provided.
  # No packages are installed automatically.

  local_library <- file.path(
    input_root,
    "local_R_library"
  )

  if (dir.exists(local_library)) {

    .libPaths(c(
      normalizePath(local_library),
      .libPaths()
    ))
  }

  required_packages <- c(
    "SummarizedExperiment",
    "ggplot2",
    "dplyr",
    "brunnermunzel",
    "survival"
  )

  available <- vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )

  if (!all(available)) {

    stop(
      "Missing R packages: ",
      paste(
        required_packages[!available],
        collapse = ", "
      )
    )
  }

  # -----------------------------------------------------
  # 3. Create a new private execution workspace.
  # -----------------------------------------------------

  stopifnot(dir.create(
    run_root,
    recursive = FALSE,
    mode = "0700"
  ))

  run_root <- normalizePath(run_root)

  # Provide access to the supplied data without copying
  # the large input matrices.
  #
  # A symbolic link is not a read-only permission
  # mechanism. The analysis scripts only read data/.

  stopifnot(file.symlink(
    file.path(input_root, "data"),
    file.path(run_root, "data")
  ))

  aggregate_dir <- file.path(
    run_root, "tables"
  )

  figure_dir <- file.path(
    run_root, "figures"
  )

  private_dir <- file.path(
    run_root, "private"
  )

  expression_private <- file.path(
    private_dir, "expression"
  )

  clinical_private <- file.path(
    private_dir, "clinical"
  )

  metastases_private <- file.path(
    private_dir, "metastases"
  )

  directories <- c(
    aggregate_dir,
    figure_dir,
    expression_private,
    clinical_private,
    metastases_private
  )

  for (directory in directories) {

    dir.create(
      directory,
      recursive = TRUE,
      showWarnings = FALSE,
      mode = "0700"
    )
  }

  # -----------------------------------------------------
  # 4. Configure all input and output paths.
  # -----------------------------------------------------

  options(

    # Data preparation.

    biomod.prepare.input_root = input_root,

    biomod.prepare.work_root = run_root,

    # Main expression analysis.

    biomod.expression.input_root = run_root,

    biomod.expression.output_dir = aggregate_dir,

    biomod.expression.private_dir = expression_private,

    # Clinical analysis.

    biomod.clinical.input_root = run_root,

    biomod.clinical.output_dir = aggregate_dir,

    biomod.clinical.private_dir = clinical_private,

    # Internal expression controls.

    biomod.controls.input_root = run_root,

    biomod.controls.output_dir = aggregate_dir,

    # Normal tissue expression.

    biomod.tissues.input_root = run_root,

    biomod.tissues.output_dir = aggregate_dir,

    # Comparison of statistical tests.

    biomod.tests.input_root = run_root,

    biomod.tests.output_dir = aggregate_dir,

    # TCGA sensitivity analysis.

    biomod.tcga.input_root = run_root,

    biomod.tcga.output_dir = aggregate_dir,

    # Primary versus metastatic acral melanoma.

    biomod.metastases.input_root = run_root,

    biomod.metastases.output_dir = aggregate_dir,

    biomod.metastases.private_dir = metastases_private,

    # Original BIOMOD survival analysis.

    biomod.survival.input_root = run_root,

    biomod.survival.output_dir = aggregate_dir,

    # Figure 1 summary statistics.

    biomod.figure1.metadata_file = file.path(
      expression_private,
      "sample_metadata.csv"
    ),

    biomod.figure1.output_dir = aggregate_dir,

    # Figure generation.

    biomod.figures.input_root = run_root,

    biomod.figures.metadata_file = file.path(
      expression_private,
      "sample_metadata.csv"
    ),

    biomod.figures.tissue_file = file.path(
      aggregate_dir,
      "normal_tissues_vs_acral.csv"
    ),

    biomod.figures.output_dir = figure_dir
  )

  # -----------------------------------------------------
  # 5. Execute the complete workflow.
  # -----------------------------------------------------

  cat("\n=== BIOMOD CSPG4 ANALYSIS ===\n")

  cat("Input root:", input_root, "\n")
  cat("Repository:", repo, "\n")
  cat("Private run directory:", run_root, "\n")

  progress <- data.frame(
    script = scripts,
    completed = FALSE,
    stringsAsFactors = FALSE
  )

  previous_wd <- getwd()

  on.exit(
    setwd(previous_wd),
    add = TRUE
  )

  # Restore the historical P-value regression check.
  historical_reference <- file.path(repo, "docs",
                                    "reference_previous_p_values.csv")
  stopifnot(file.exists(historical_reference))
  historical_dir <- file.path(run_root, "results_v2/tables")
  dir.create(historical_dir, recursive = TRUE,
             showWarnings = FALSE, mode = "0700")
  historical_target <- file.path(historical_dir,
                                 "contrastes_finales_CSPG4.csv")
  stopifnot(file.copy(historical_reference, historical_target,
                     overwrite = FALSE))
  Sys.chmod(historical_target, mode = "0600")

  for (i in seq_along(script_paths)) {

    cat(
      "\n====================================================\n",
      "STEP ", i, " OF ", length(scripts),
      ": ", scripts[i],
      "\n====================================================\n",
      sep = ""
    )

    setwd(repo)

    # Each script has an independent evaluation
    # environment and receives explicit path options.
    #
    # Statistical scripts set their verified RNG settings.
    # Script 14 separately sets its graphics RNG settings.

    source(
      script_paths[i],
      local = new.env(parent = globalenv())
    )

    progress$completed[i] <- TRUE

    write.csv(
      progress,
      file.path(private_dir, "workflow_progress.csv"),
      row.names = FALSE
    )

    invisible(gc())

    cat("\nCompleted:", scripts[i], "\n")
  }

  # -----------------------------------------------------
  # 6. Verify the expected final output inventory.
  # -----------------------------------------------------

  expected_tables <- c(

    "acral_metastasis_overlap_sensitivity.csv",

    "acral_primary_vs_lymph_node_metastases.csv",

    "acral_survival_biomod.csv",

    "acral65_sensitivity.csv",

    "analysis_checks.csv",

    "clinical_association_tests.csv",

    "clinical_cluster_pairwise.csv",

    "clinical_descriptive.csv",

    "CSPG4_vs_deconvolution.csv",

    "CSPG4_vs_purity.csv",

    "dunn_sensitivity.csv",

    "expression_controls_contrasts.csv",

    "expression_controls_summary.csv",

    "expression_descriptive.csv",

    "figure1_summary_statistics.csv",

    "main_contrasts.csv",

    "normal_tissue_expression_summary.csv",

    "normal_tissues_vs_acral.csv",

    "statistical_test_comparison.csv",

    "TCGA_sensitivity_analysis.csv",

    "vascular_marker_correlations.csv"
  )

  expected_figures <- c(

    "CSPG4_expression_comparison.png",

    "CSPG4_normal_tissues_vs_acral.png",

    "CSPG4_acral44_clinical.png"
  )

  actual_tables <- list.files(
    aggregate_dir,
    pattern = "\\.csv$"
  )

  actual_figures <- list.files(
    figure_dir,
    pattern = "\\.png$"
  )

  stopifnot(
    all(progress$completed),
    setequal(actual_tables, expected_tables),
    setequal(actual_figures, expected_figures)
  )

  main_checks <- read.csv(file.path(
    aggregate_dir,
    "analysis_checks.csv"
  ))

  stopifnot(
    "ok" %in% names(main_checks),
    nrow(main_checks) == 17L,
    all(main_checks$ok)
  )

  # Record the actual execution environment privately.

  writeLines(
    capture.output(sessionInfo()),
    file.path(private_dir, "sessionInfo.txt")
  )

  write.csv(
    data.frame(
      file = c(expected_tables, expected_figures),
      type = c(
        rep("Aggregate table", length(expected_tables)),
        rep("Figure", length(expected_figures))
      )
    ),
    file.path(private_dir, "output_inventory.csv"),
    row.names = FALSE
  )

  cat("\n=== COMPLETE WORKFLOW STATUS ===\n")

  print(progress, row.names = FALSE)

  cat(
    "\nScripts completed:",
    sum(progress$completed),
    "/", nrow(progress),
    "\n"
  )

  cat(
    "Aggregate CSV tables:",
    length(actual_tables),
    "\n"
  )

  cat(
    "PNG figures:",
    length(actual_figures),
    "\n"
  )

  cat(
    "Main analysis checks passed:",
    sum(main_checks$ok),
    "/", nrow(main_checks),
    "\n"
  )

  cat(
    "\nSUCCESS: Complete BIOMOD workflow executed.\n"
  )

  cat(
    "\nAll generated files remain in:\n",
    run_root,
    "\n"
  )

  cat(
    "\nReview data-sharing authorization before ",
    "distributing tables or figures.\n",
    sep = ""
  )

}

main()
