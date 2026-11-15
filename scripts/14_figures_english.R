# 14_figures_english.R
# BIOMOD 2026 | Pandora's Box | Objective 1
# Regenerate the three CSPG4 figures in English from the existing analysis.
# Configure input and output options before sourcing this script.
# No statistical tests, cohort selections or existing CSV files are modified.
# Write three English PNG files to the configured output directory.
# Individual-level inputs must remain local and must not be uploaded to GitHub.

local({
  input_root <- getOption("biomod.figures.input_root", ".")

  metadata_file <- getOption(
    "biomod.figures.metadata_file",
    file.path(
      getOption("biomod.expression.private_dir",
                file.path(tempdir(), "biomod_private")),
      "sample_metadata.csv"
    )
  )

  tissue_file <- getOption(
    "biomod.figures.tissue_file",
    file.path(
      getOption("biomod.expression.output_dir", "tables"),
      "normal_tissues_vs_acral.csv"
    )
  )

  figure_dir <- getOption("biomod.figures.output_dir")

  if (is.null(figure_dir)) {
    stop("Set biomod.figures.output_dir explicitly.")
  }

  # Figures 1 and 3 contain individual expression points.
  # Review data-sharing authorization before publication.

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 first: install.packages('ggplot2')")
  }
  suppressPackageStartupMessages(library(ggplot2))
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(20261022)

  # Existing outputs and the same local clinical inputs as scripts 05 and 06.
  input <- c(
    expression = metadata_file,
    tissues = tissue_file,
    acral = file.path(input_root, "results_v2", "tables",
                      "nivelA_CSPG4_TPMcrudo_86muestras.csv"),
    clinical = file.path(input_root, "data", "raw",
                         "acralMX_metadata_cc_update.csv"),
    clusters = file.path(input_root, "results_v2", "tables",
                         "nivelA_CSPG4_44muestras.csv")
  )

  missing_files <- input[!file.exists(input)]
  if (length(missing_files)) {
    stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
  }
  require_columns <- function(x, columns, label) {
    missing <- setdiff(columns, names(x))
    if (length(missing)) stop(label, ": missing column(s): ", paste(missing, collapse = ", "))
  }

  m <- read.csv(input[["expression"]], check.names = FALSE)
  tissues <- read.csv(input[["tissues"]], check.names = FALSE)
  a <- read.csv(input[["acral"]], check.names = FALSE)
  clinical <- read.csv(input[["clinical"]], check.names = FALSE)
  clusters <- read.csv(input[["clusters"]], check.names = FALSE)
  require_columns(m, c("group", "log2TPM1"), "Expression table")
  require_columns(tissues, c("tissue", "median_TPM", "log2TPM1"), "Tissue table")
  require_columns(a, c("en_44", "sample_id", "CSPG4_log2"), "Acral table")
  require_columns(clinical, c("sample_id", "stage", "ulceration", "sex"), "Clinical table")
  require_columns(clusters, "sample_id", "Cluster table")

  # ---- Validate the existing cohorts (no new filtering) -------------------
  group_codes <- c("G1a", "G1b", "G2", "G3")
  m$group <- factor(sub(" .*", "", as.character(m$group)), levels = group_codes)
  stopifnot(!anyNA(m$group),
            identical(as.integer(table(m$group)), c(604L, 701L, 103L, 44L)),
            all(is.finite(m$log2TPM1)))
  values <- split(m$log2TPM1, m$group)
  observed_medians <- vapply(values, median, numeric(1))
  stopifnot(all(abs(observed_medians - c(3.5759, 3.8298, 7.6731, 5.5107)) < 1e-3))
  acral_median <- median(values[["G3"]])

  stopifnot(nrow(tissues) == 54L, !anyDuplicated(tissues$tissue),
            all(is.finite(tissues$median_TPM)), all(is.finite(tissues$log2TPM1)),
            !anyNA(tissues$tissue))
  stopifnot(sum(tissues$median_TPM > median(2^values[["G3"]] - 1)) == 10L)

  # Clinical reconstruction follows the original script 06.
  # No patient identifiers are included in the plotting data or output.
  a44 <- a[a$en_44, c("sample_id", "CSPG4_log2")]
  stopifnot(nrow(a44) == 44L, !anyDuplicated(a44$sample_id))
  d <- merge(a44, clinical, by = "sample_id")
  stopifnot(nrow(d) == 44L, !anyDuplicated(d$sample_id),
            all(is.finite(d$CSPG4_log2)),
            isTRUE(all.equal(sort(d$CSPG4_log2), sort(values[["G3"]]), tolerance = 1e-8)))
  cluster_column <- grep("cluster", names(clusters), ignore.case = TRUE, value = TRUE)
  if (!length(cluster_column)) stop("Cluster column not found in the 44-sample table.")
  stopifnot(!anyDuplicated(clusters$sample_id))
  d$cluster <- factor(clusters[[cluster_column[1]]][match(d$sample_id, clusters$sample_id)])
  stopifnot(!anyNA(d$cluster), identical(sort(as.integer(table(d$cluster))), c(14L, 14L, 16L)))
  stopifnot(all(na.omit(as.character(d$ulceration)) %in% c("0", "1")),
            all(na.omit(as.character(d$sex)) %in% c("F", "M")),
            sum(!is.na(d$ulceration)) == 43L,
            sum(!is.na(d$stage)) == 44L,
            sum(!is.na(d$sex)) == 44L)

  # ---- Figure 1: comparison of expression across four groups --------------
  # The original violin, boxplot, jitter and default group palette are retained.
  labels <- setNames(c(
    "G1a\nNon-sun-exposed skin\n(GTEx; n = 604)",
    "G1b\nSun-exposed skin\n(GTEx; n = 701)",
    "G2\nPrimary melanoma\n(TCGA-SKCM; n = 103)",
    "G3\nPrimary acral melanoma\n(Mexico; n = 44)"
  ), group_codes)
  p1 <- ggplot(m, aes(group, log2TPM1, fill = group)) +
    geom_violin(alpha = .35, colour = NA, scale = "width") +
    geom_boxplot(width = .18, outlier.shape = NA, fill = "white") +
    geom_jitter(width = .08, size = .6, alpha = .35) +
    scale_x_discrete(labels = labels) +
    labs(x = NULL, y = "CSPG4 expression [log2(TPM + 1)]",
         title = "CSPG4 expression in non-acral skin and primary melanoma",
         subtitle = "Exploratory cross-study comparison; no additional batch-effect correction") +
    theme_bw(base_size = 11) + theme(legend.position = "none")

  # ---- Figure 2: expression across GTEx tissue categories ------------------
  # The 54 original GTEx names are retained; two categories are cell lines.
  tissues$tissue <- factor(tissues$tissue,
                           levels = rev(tissues$tissue[order(-tissues$median_TPM)]))
  tissues$type <- ifelse(grepl("^Skin", tissues$tissue), "Skin", "Other tissues")
  p2 <- ggplot(tissues, aes(tissue, log2TPM1, fill = type)) +
    geom_col() +
    geom_hline(yintercept = acral_median, colour = "firebrick", linetype = "dashed") +
    scale_fill_manual(values = c("Skin" = "steelblue", "Other tissues" = "grey65")) +
    coord_flip() +
    labs(x = NULL, y = "Median CSPG4 expression [log2(TPM + 1)]", fill = NULL,
         title = "CSPG4 expression across 54 GTEx v8 tissue categories",
         subtitle = sprintf(paste0("Red dashed line: primary acral melanoma median (%.2f)\n",
                                   "Bulk-tissue medians; descriptive comparison"), acral_median)) +
    theme_bw(base_size = 9)

  # ---- Figure 3: clinical and molecular characteristics --------------------
  # Preserve the original panel order (Cluster, Stage, Sex, Ulceration).
  # Translate display labels only; keep the underlying categories unchanged.
  d$ulc <- factor(d$ulceration, levels = c(0, 1),
                  labels = c("without_ulceration", "with_ulceration"))
  long <- rbind(
    data.frame(panel = "Cluster", level = as.character(d$cluster), y = d$CSPG4_log2),
    data.frame(panel = "Clinical stage", level = as.character(d$stage), y = d$CSPG4_log2),
    data.frame(panel = "Sex", level = as.character(d$sex), y = d$CSPG4_log2),
    data.frame(panel = "Ulceration", level = as.character(d$ulc), y = d$CSPG4_log2)
  )
  long <- long[!is.na(long$level), ]
  long$panel <- factor(long$panel, levels = c("Cluster", "Clinical stage", "Sex", "Ulceration"))
  stopifnot(identical(as.integer(table(long$panel)), c(44L, 44L, 44L, 43L)))
  display_clinical <- function(x) {
    map <- c(F = "Female", M = "Male", with_ulceration = "Ulcerated",
             without_ulceration = "Non-ulcerated")
    idx <- match(x, names(map))
    use <- !is.na(idx)
    x[use] <- unname(map[idx[use]])
    x
  }
  p3 <- ggplot(long, aes(level, y)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90") +
    geom_jitter(width = .15, size = 1.6, alpha = .7) +
    facet_wrap(~panel, scales = "free_x", nrow = 1) +
    scale_x_discrete(labels = display_clinical) +
    labs(x = NULL, y = "CSPG4 expression [log2(TPM + 1)]",
         title = "CSPG4 expression by clinical and molecular characteristics",
         subtitle = "44 primary acral melanomas; exploratory, small sample, no confounder adjustment") +
    theme_bw(base_size = 11)

  # ---- Save only new English figures --------------------------------------
  out <- figure_dir
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(out, c("CSPG4_expression_comparison.png",
                            "CSPG4_normal_tissues_vs_acral.png",
                            "CSPG4_acral44_clinical.png"))
  if (any(file.exists(files))) {
    stop("Figure files already exist. Use a new output directory.")
  }

  ggsave(files[1], p1, width = 9.5, height = 5.8, dpi = 300)
  ggsave(files[2], p2, width = 7.5, height = 9, dpi = 300)
  ggsave(files[3], p3, width = 11, height = 4, dpi = 200)

  # Print aggregate checks only. These do not overwrite any statistical output.
  cat("\nEnglish figures saved:\n", paste(files, collapse = "\n"), "\n", sep = "")
  cat("\nSample counts: ", paste(as.integer(table(m$group)), collapse = "/"), "\n", sep = "")
  cat("Acral melanoma median: ", round(acral_median, 3), " log2(TPM + 1)\n", sep = "")
  cat("Original figures and analysis tables were not overwritten.\n")
})
