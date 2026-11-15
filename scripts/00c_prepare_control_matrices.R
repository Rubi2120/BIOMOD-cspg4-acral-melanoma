# BIOMOD 2026 - Pandora's Box
# Prepare complete TPM matrices for expression controls.
#
# G1a: GTEx non-sun-exposed skin.
# G1b: GTEx sun-exposed skin.
# G2: TCGA primary melanoma.
# G3: Mexican primary acral melanoma.
#
# All matrices contain original TPM values.
# Ensembl version suffixes are removed and the first
# occurrence of duplicate gene identifiers is retained,
# following the original BIOMOD implementation.
#
# This script does not perform batch correction or
# normalization across studies.
#
# The resulting RDS contains individual-level expression
# data and must remain private.

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

  input <- function(path) {
    file.path(input_root, path)
  }

  table_dir <- file.path(
    work_root, "results_v2/tables"
  )

  output_dir <- file.path(
    work_root, "results_v2/intermedios"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE,
    mode = "0700"
  )

  output_file <- file.path(
    output_dir, "ctl_tpm.rds"
  )

  if (file.exists(output_file)) {
    stop("Control matrix file already exists. Use a new private workspace.")
  }

  strip_version <- function(x) {
    sub("\\..*$", "", x)
  }

  remove_duplicates <- function(m) {
    m[!duplicated(rownames(m)), , drop = FALSE]
  }

  genes <- c(
    CSPG4 = "ENSG00000173546",
    MLANA = "ENSG00000120215",
    PMEL  = "ENSG00000185664",
    ACTB  = "ENSG00000075624",
    GAPDH = "ENSG00000111640",
    TBP   = "ENSG00000112592"
  )

  # Read the sample-selection tables produced by script 00.

  a <- read.csv(file.path(
    table_dir,
    "nivelA_CSPG4_TPMcrudo_86muestras.csv"
  ))

  b <- read.csv(file.path(
    table_dir,
    "nivelB_TCGA_CSPG4_por_muestra.csv"
  ))

  p <- read.csv(file.path(
    table_dir,
    "nivelC_GTEx_CSPG4_piel_muestras.csv"
  ))

  # -------------------------------------------------------
  # G3. Mexican acral melanoma: 44 primary tumors
  # -------------------------------------------------------

  t86 <- read.csv(
    input("data/raw/acral_MX_TPMs_qcpassed.tsv"),
    check.names = FALSE
  )

  ids44 <- a$sample_id[a$en_44]

  stopifnot(
    length(ids44) == 44,
    all(ids44 %in% colnames(t86))
  )

  G3 <- as.matrix(t86[, ids44])

  storage.mode(G3) <- "double"

  rownames(G3) <- strip_version(
    t86$ENSEMBL_GENE_ID
  )

  G3 <- remove_duplicates(G3)

  rm(t86)
  invisible(gc())

  # -------------------------------------------------------
  # G2. TCGA-SKCM primary tumors
  # -------------------------------------------------------

  se <- readRDS(
    input("data/raw/GDCdata/tcga_skcm_se_preparado.rds")
  )

  gid <- if ("gene_id" %in% names(rowData(se))) {
    rowData(se)$gene_id
  } else {
    rownames(se)
  }

  ids2 <- b$barcode[
    b$sample_type == "Primary Tumor"
  ]

  stopifnot(length(ids2) == 103)

  G2 <- as.matrix(
    assay(se, "tpm_unstrand")[, ids2]
  )

  rownames(G2) <- strip_version(gid)

  G2 <- remove_duplicates(G2)

  rm(se)
  invisible(gc())

  # -------------------------------------------------------
  # G1a and G1b. GTEx skin samples
  # -------------------------------------------------------

  non_sun <- grepl(
    "Not Sun Exposed", p$SMTSD
  )

  id1a <- p$SAMPID[non_sun]

  id1b <- p$SAMPID[!non_sun]

  tpm_file <- input(
    paste0(
      "data/raw/",
      "GTEx_Analysis_2017-06-05_v8_",
      "RNASeQCv1.1.9_gene_tpm.gct.gz"
    )
  )

  # Stream the original GCT file in 500-gene blocks.
  # Preserve the original sample and gene ordering.

  read_skin_matrix <- function(path, sample_ids) {

    con <- gzfile(path, "r")

    on.exit(close(con))

    invisible(readLines(con, 2))

    header <- strsplit(
      readLines(con, 1),
      "\t",
      fixed = TRUE
    )[[1]]

    selected <- match(
      sample_ids, header
    )

    stopifnot(!anyNA(selected))

    matrices <- list()
    gene_ids <- list()
    descriptions <- list()

    k <- 0

    message("Reading the complete GTEx skin matrix.")

    repeat {

      block <- readLines(con, 500)

      if (!length(block)) break

      fields <- strsplit(
        block, "\t", fixed = TRUE
      )

      k <- k + 1

      gene_ids[[k]] <- vapply(
        fields, `[`, "", 1
      )

      descriptions[[k]] <- vapply(
        fields, `[`, "", 2
      )

      matrices[[k]] <- t(vapply(
        fields,
        function(z) as.numeric(z[selected]),
        numeric(length(selected))
      ))

      if (k %% 10 == 0) {
        message("  Genes read: ", k * 500)
      }
    }

    M <- do.call(rbind, matrices)

    rownames(M) <- strip_version(
      unlist(gene_ids)
    )

    colnames(M) <- header[selected]

    desc <- unlist(descriptions)

    duplicate <- duplicated(rownames(M))

    list(
      matrix = M[!duplicate, , drop = FALSE],
      description = desc[!duplicate]
    )
  }

  skin <- read_skin_matrix(
    tpm_file, c(id1a, id1b)
  )

  M <- skin$matrix

  desc <- skin$description

  G1a <- M[, id1a, drop = FALSE]

  G1b <- M[, id1b, drop = FALSE]

  tpm <- list(
    G1a = G1a,
    G1b = G1b,
    G2 = G2,
    G3 = G3
  )

  # -------------------------------------------------------
  # Verification against the prepared expression tables
  # -------------------------------------------------------

  n <- vapply(tpm, ncol, numeric(1))

  stopifnot(identical(
    unname(n),
    c(604, 701, 103, 44)
  ))

  symbols <- desc[
    match(genes, rownames(M))
  ]

  stopifnot(identical(
    unname(symbols),
    names(genes)
  ))

  for (group in names(tpm)) {

    missing_genes <- setdiff(
      genes, rownames(tpm[[group]])
    )

    if (length(missing_genes)) {
      stop("Missing control genes in ", group)
    }
  }

  cs <- genes[["CSPG4"]]

  checks <- c(

    G1 = max(abs(
      log2(M[cs, p$SAMPID] + 1) -
      p$CSPG4_log2TPM1
    )),

    G2 = max(abs(
      log2(
        G2[
          cs,
          b$barcode[b$sample_type == "Primary Tumor"]
        ] + 1
      ) -
      b$CSPG4_log2TPM1[
        b$sample_type == "Primary Tumor"
      ]
    )),

    G3 = max(abs(
      log2(G3[cs, ids44] + 1) -
      a$CSPG4_log2[a$en_44]
    ))
  )

  cat(
    "\nMaximum CSPG4 differences against expression tables:\n"
  )

  print(checks)

  stopifnot(all(checks < 1e-6))

  common_genes <- Reduce(
    intersect,
    lapply(tpm, rownames)
  )

  stopifnot(
    length(common_genes) == 55598
  )

  # GTEx tissue median matrix, used by script 08.

  med <- read.delim(
    gzfile(input(
      paste0(
        "data/raw/",
        "GTEx_Analysis_2017-06-05_v8_",
        "RNASeQCv1.1.9_gene_median_tpm.gct.gz"
      )
    )),
    skip = 2,
    check.names = FALSE
  )

  stopifnot(ncol(med) == 56)

  # Preserve the original intermediate RDS structure.

  saveRDS(
    list(
      tpm = tpm,
      genes = genes,
      med = med
    ),
    output_file
  )

  cat("\nMatrix dimensions:\n")

  print(vapply(tpm, dim, integer(2)))

  cat(
    "\nGenes shared by all four groups:",
    length(common_genes),
    "\n"
  )

  cat(
    "\nPrivate control matrices saved to:\n",
    output_file,
    "\n"
  )

})
