# BIOMOD 2026 - Pandora's Box
# Prepare CSPG4 expression tables from original inputs.
#
# All expression values are obtained from previously
# quantified TPM matrices. No RNA sequencing or read-level
# processing is performed here.
#
# The 44-sample acral expression matrix is used to identify
# the main analysis subset. Expression values for the main
# analysis are taken from the original 86-sample TPM matrix.
#
# Original intermediate filenames and column names are
# retained for compatibility with scripts 05 through 12.
#
# All intermediate outputs must remain private.

local({

  suppressPackageStartupMessages(
    library(SummarizedExperiment)
  )

  input_root <- getOption(
    "biomod.prepare.input_root", "."
  )

  work_root <- getOption(
    "biomod.prepare.work_root"
  )

  if (is.null(work_root)) {
    stop("Set biomod.prepare.work_root to a private directory.")
  }

  output_dir <- file.path(
    work_root, "results_v2/tables"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE,
    mode = "0700"
  )

  input <- function(path) {
    file.path(input_root, path)
  }

  output <- function(path) {
    file.path(output_dir, path)
  }

  filenames <- c(
    "nivelA_CSPG4_TPMcrudo_86muestras.csv",
    "nivelB_TCGA_CSPG4_por_muestra.csv",
    "GTEx_CSPG4_mediana_54tejidos.csv",
    "nivelC_GTEx_CSPG4_piel_muestras.csv"
  )

  if (any(file.exists(output(filenames)))) {
    stop("Intermediate tables already exist. Use a new private workspace.")
  }

  strip_version <- function(x) {
    sub("\\..*$", "", x)
  }

  ENSG <- "ENSG00000173546"

  # -------------------------------------------------------
  # A. Mexican acral melanoma: original TPM values
  # -------------------------------------------------------

  # Despite its .tsv extension, this file is comma-separated.

  tpm86 <- read.csv(
    input("data/raw/acral_MX_TPMs_qcpassed.tsv"),
    check.names = FALSE
  )

  meta86 <- read.csv(
    input("data/raw/acral_MX_metadata_.csv"),
    check.names = FALSE
  )

  expr44 <- readRDS(
    input("data/raw/acralMX_expression_data_update44s.rds")
  )

  i <- which(
    strip_version(tpm86$ENSEMBL_GENE_ID) == ENSG
  )

  stopifnot(length(i) == 1)

  samples <- colnames(tpm86)[-(1:2)]

  raw_cspg4 <- as.numeric(
    tpm86[i, samples]
  )

  ids44 <- strip_version(
    rownames(expr44)
  )

  # Preserve the processed expression values for
  # provenance and comparison only.

  processed_cspg4 <- as.numeric(
    expr44[match(ENSG, ids44), ][
      match(samples, colnames(expr44))
    ]
  )

  acral <- data.frame(
    sample_id = samples,
    sample_type = meta86$sample_type[
      match(samples, meta86$sample_id)
    ],
    en_44 = samples %in% colnames(expr44),
    CSPG4_TPM = raw_cspg4,
    CSPG4_log2 = log2(raw_cspg4 + 1),
    CSPG4_matriz44_recortada = processed_cspg4,
    prefijo = substr(samples, 1, 7)
  )

  stopifnot(
    nrow(acral) == 86,
    sum(acral$en_44) == 44,
    !anyNA(acral$sample_type),
    !anyDuplicated(acral$sample_id)
  )

  write.csv(
    acral,
    output(filenames[1]),
    row.names = FALSE
  )

  rm(tpm86, expr44)
  invisible(gc())

  # -------------------------------------------------------
  # B. TCGA-SKCM: GDC unstranded TPM
  # -------------------------------------------------------

  se <- readRDS(
    input("data/raw/GDCdata/tcga_skcm_se_preparado.rds")
  )

  ids <- if ("gene_id" %in% names(rowData(se))) {
    rowData(se)$gene_id
  } else {
    rownames(se)
  }

  i <- which(
    strip_version(ids) == ENSG
  )

  stopifnot(length(i) == 1)

  cd <- as.data.frame(colData(se))

  tpm <- as.numeric(
    assay(se, "tpm_unstrand")[i, ]
  )

  tcga <- data.frame(
    barcode = colnames(se),
    paciente = substr(colnames(se), 1, 12),
    sample_type = cd$sample_type,
    CSPG4_TPM = tpm,
    CSPG4_log2TPM1 = log2(tpm + 1)
  )

  stopifnot(
    nrow(tcga) == 471,
    sum(tcga$sample_type == "Primary Tumor") == 103,
    !anyDuplicated(tcga$barcode)
  )

  write.csv(
    tcga,
    output(filenames[2]),
    row.names = FALSE
  )

  rm(se)
  invisible(gc())

  # -------------------------------------------------------
  # GTEx: tissue median expression
  # -------------------------------------------------------

  median_file <- input(
    paste0(
      "data/raw/",
      "GTEx_Analysis_2017-06-05_v8_",
      "RNASeQCv1.1.9_gene_median_tpm.gct.gz"
    )
  )

  med <- read.delim(
    gzfile(median_file),
    skip = 2,
    check.names = FALSE
  )

  i <- which(
    strip_version(med$Name) == ENSG
  )

  stopifnot(length(i) == 1)

  tissues <- data.frame(
    tejido = colnames(med)[-(1:2)],
    mediana_TPM = as.numeric(med[i, -(1:2)])
  )

  tissues$log2TPM1 <- log2(
    tissues$mediana_TPM + 1
  )

  tissues <- tissues[
    order(-tissues$mediana_TPM),
  ]

  stopifnot(nrow(tissues) == 54)

  write.csv(
    tissues,
    output(filenames[3]),
    row.names = FALSE
  )

  # -------------------------------------------------------
  # C. GTEx: individual skin samples
  # -------------------------------------------------------

  tpm_file <- input(
    paste0(
      "data/raw/",
      "GTEx_Analysis_2017-06-05_v8_",
      "RNASeQCv1.1.9_gene_tpm.gct.gz"
    )
  )

  # Read the compressed GCT file in blocks and retain
  # only the CSPG4 expression row.

  read_gene <- function(path, gene_id) {

    con <- gzfile(path, "r")

    on.exit(close(con))

    invisible(readLines(con, 2))

    header <- strsplit(
      readLines(con, 1), "\t"
    )[[1]]

    gene_row <- NULL

    repeat {

      block <- readLines(con, 2000)

      if (!length(block)) break

      hits <- grep(
        paste0("^", gene_id, "\\."),
        block
      )

      if (length(hits)) {

        gene_row <- strsplit(
          block[hits[1]], "\t"
        )[[1]]

        break
      }
    }

    stopifnot(!is.null(gene_row))

    setNames(
      as.numeric(gene_row[-(1:2)]),
      header[-(1:2)]
    )
  }

  gtex_tpm <- read_gene(
    tpm_file, ENSG
  )

  skin_meta <- read.csv(
    input("data/GTEx_skin_samples_v8.csv")
  )

  skin_meta <- skin_meta[
    skin_meta$SAMPID %in% names(gtex_tpm),
  ]

  # Preserve the exact five-column schema of the
  # original BIOMOD analysis table.
  #
  # The original preparation script generated an extra
  # exposure-site column. It is omitted here because
  # exposure status is already encoded in SMTSD.

  skin <- data.frame(
    SAMPID = skin_meta$SAMPID,
    SMTSD = skin_meta$SMTSD,
    CSPG4_TPM = unname(
      gtex_tpm[skin_meta$SAMPID]
    ),
    donante = sub(
      "^(GTEX-[^-]+)-.*$",
      "\\1",
      skin_meta$SAMPID
    )
  )

  skin$CSPG4_log2TPM1 <- log2(
    skin$CSPG4_TPM + 1
  )

  non_sun <- grepl(
    "Not Sun Exposed",
    skin$SMTSD
  )

  stopifnot(
    nrow(skin) == 1305,
    sum(non_sun) == 604,
    sum(!non_sun) == 701,
    !anyNA(skin$CSPG4_TPM),
    !anyDuplicated(skin$SAMPID)
  )

  write.csv(
    skin,
    output(filenames[4]),
    row.names = FALSE
  )

  # Report only aggregate sample counts.

  cat("\nPrepared expression tables:\n")

  print(data.frame(
    table = filenames,
    rows = c(
      nrow(acral),
      nrow(tcga),
      nrow(tissues),
      nrow(skin)
    )
  ), row.names = FALSE)

  cat(
    "\nPrivate intermediate tables saved to:\n",
    output_dir, "\n"
  )

})
