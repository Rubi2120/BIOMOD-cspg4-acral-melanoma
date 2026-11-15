# BIOMOD 2026 - Pandora's Box
# Main analysis of CSPG4 expression in non-acral skin,
# primary TCGA melanoma and Mexican primary acral melanoma.
#
# This script reproduces the statistical analysis reported
# at BIOMOD. The contrasts were defined during data exploration
# and should be interpreted as exploratory.
#
# Individual-level inputs and outputs must remain private.
# The final English figures are generated separately by script 14.

input_root <- getOption("biomod.expression.input_root", ".")
output_dir <- getOption("biomod.expression.output_dir", "tables")
private_dir <- getOption(
  "biomod.expression.private_dir",
  file.path(tempdir(), "biomod_private")
)

# Restore the original BIOMOD random number configuration.
RNGkind(
  kind = "L'Ecuyer-CMRG",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

set.seed(20261022)
B <- 2000

# Preserve the original package initialization sequence.
for (p in c("ggplot2", "brunnermunzel")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Required package: ", p)
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)

checks <- list()

add_check <- function(test, detail, ok) {
  checks[[length(checks) + 1]] <<- data.frame(
    test = test,
    detail = as.character(detail),
    ok = isTRUE(ok)
  )
}

# 1. Input data -------------------------------------------------

input_paths <- c(
  acral = "results_v2/tables/nivelA_CSPG4_TPMcrudo_86muestras.csv",
  tcga = "results_v2/tables/nivelB_TCGA_CSPG4_por_muestra.csv",
  gtex = "results_v2/tables/nivelC_GTEx_CSPG4_piel_muestras.csv",
  tissues = "results_v2/tables/GTEx_CSPG4_mediana_54tejidos.csv"
)

input_paths <- setNames(
  file.path(input_root, input_paths),
  names(input_paths)
)

missing_inputs <- input_paths[!file.exists(input_paths)]

if (length(missing_inputs)) {
  stop("Missing inputs:\n", paste(missing_inputs, collapse = "\n"))
}

a <- read.csv(input_paths[["acral"]])
tc <- read.csv(input_paths[["tcga"]])
gt <- read.csv(input_paths[["gtex"]])
tissue <- read.csv(input_paths[["tissues"]])

# 2. Integrated sample-level dataset ---------------------------

groups <- c(
  "G1a GTEx non-sun-exposed skin",
  "G1b GTEx sun-exposed skin",
  "G2 TCGA primary cutaneous melanoma",
  "G3 Mexican primary acral melanoma"
)

gt$site <- ifelse(
  grepl("Not Sun", gt$SMTSD),
  "non_sun_exposed",
  "sun_exposed"
)

stopifnot(!anyDuplicated(paste(gt$donante, gt$site)))

tc1 <- tc[tc$sample_type == "Primary Tumor", ]
a44 <- a[a$en_44, ]
a65 <- a[a$sample_type == "primary", ]

m <- rbind(
  data.frame(
    sample = gt$SAMPID[gt$site == "non_sun_exposed"],
    group = groups[1],
    source = "GTEx v8",
    platform = "RNA-seq, RNA-SeQC 1.1.9",
    unit = "TPM",
    log2TPM1 = gt$CSPG4_log2TPM1[
      gt$site == "non_sun_exposed"
    ]
  ),
  data.frame(
    sample = gt$SAMPID[gt$site == "sun_exposed"],
    group = groups[2],
    source = "GTEx v8",
    platform = "RNA-seq, RNA-SeQC 1.1.9",
    unit = "TPM",
    log2TPM1 = gt$CSPG4_log2TPM1[
      gt$site == "sun_exposed"
    ]
  ),
  data.frame(
    sample = tc1$barcode,
    group = groups[3],
    source = "TCGA-SKCM (GDC)",
    platform = "RNA-seq, GDC (tpm_unstrand)",
    unit = "TPM",
    log2TPM1 = tc1$CSPG4_log2TPM1
  ),
  data.frame(
    sample = a44$sample_id,
    group = groups[4],
    source = "Mexican acral melanoma cohort",
    platform = "RNA-seq, exome capture",
    unit = "TPM",
    log2TPM1 = a44$CSPG4_log2
  )
)

m$group <- factor(m$group, levels = groups)

# Private output: never distribute individual-level data.
write.csv(
  m,
  file.path(private_dir, "sample_metadata.csv"),
  row.names = FALSE
)

v <- split(m$log2TPM1, m$group)

add_check(
  "Sample counts: 604/701/103/44",
  paste(lengths(v), collapse = "/"),
  identical(as.integer(lengths(v)), c(604L, 701L, 103L, 44L))
)

add_check(
  "Unique sample identifiers",
  sum(duplicated(m$sample)),
  !anyDuplicated(m$sample)
)

add_check(
  "Finite expression values",
  sum(!is.finite(m$log2TPM1)),
  all(is.finite(m$log2TPM1))
)

add_check(
  "TCGA primary tumors: unique patients",
  sum(duplicated(tc1$paciente)),
  !anyDuplicated(tc1$paciente)
)

add_check(
  "Acral: 44 primary tumors among 86 samples",
  paste(nrow(a44), unique(a44$sample_type)),
  nrow(a44) == 44 && all(a44$sample_type == "primary")
)

add_check(
  "Acral: 65 primary tumors",
  nrow(a65),
  nrow(a65) == 65
)

# 3. Descriptive statistics ------------------------------------

desc <- do.call(rbind, lapply(names(v), function(g) {
  data.frame(
    group = g,
    n = length(v[[g]]),
    mean = mean(v[[g]]),
    median = median(v[[g]]),
    Q1 = quantile(v[[g]], .25, names = FALSE),
    Q3 = quantile(v[[g]], .75, names = FALSE)
  )
}))

write.csv(
  desc,
  file.path(output_dir, "expression_descriptive.csv"),
  row.names = FALSE
)

add_check(
  "Expected group medians",
  paste(round(desc$median, 4), collapse = "/"),
  all(abs(
    desc$median - c(3.5759, 3.8298, 7.6731, 5.5107)
  ) < 1e-3)
)

# 4. Main statistical contrasts --------------------------------

# Cliff's delta: P(x > y) - P(x < y), accounting for ties.
cliff_delta <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  r <- rank(c(x, y))

  2 * ((sum(r[1:n1]) - n1 * (n1 + 1) / 2) /
         (n1 * n2)) - 1
}

# Bootstrap confidence intervals.
bootstrap_ci <- function(x, y, B) {
  d <- replicate(B, {
    xb <- sample(x, replace = TRUE)
    yb <- sample(y, replace = TRUE)

    c(
      cliff_delta(xb, yb),
      median(xb) - median(yb)
    )
  })

  apply(d, 1, quantile, probs = c(.025, .975))
}

compare_groups <- function(g1, g2) {

  x <- v[[g1]]
  y <- v[[g2]]

  bm <- brunnermunzel::brunnermunzel.test(x, y)
  w <- suppressWarnings(wilcox.test(x, y))
  ci <- bootstrap_ci(x, y, B)

  data.frame(
    contrast = paste(g1, "vs", g2),
    method = "Brunner-Munzel",
    n1 = length(x),
    n2 = length(y),
    HL = NA_real_,
    HL_CI_lower = NA_real_,
    HL_CI_upper = NA_real_,
    cliff_delta = cliff_delta(x, y),
    delta_CI_lower = ci[1, 1],
    delta_CI_upper = ci[2, 1],
    median_difference = median(x) - median(y),
    median_diff_CI_lower = ci[1, 2],
    median_diff_CI_upper = ci[2, 2],
    p = bm$p.value,
    p_wilcoxon = w$p.value
  )
}

# 4a. Paired comparison of GTEx skin samples.

w <- reshape(
  gt[, c("donante", "site", "CSPG4_log2TPM1")],
  idvar = "donante",
  timevar = "site",
  direction = "wide"
)

w <- w[complete.cases(w), ]

wp <- suppressWarnings(wilcox.test(
  w$CSPG4_log2TPM1.sun_exposed,
  w$CSPG4_log2TPM1.non_sun_exposed,
  paired = TRUE,
  conf.int = TRUE
))

paired_result <- data.frame(
  contrast = paste(groups[2], "vs", groups[1], "(paired)"),
  method = "Paired Wilcoxon",
  n1 = nrow(w),
  n2 = nrow(w),
  HL = unname(wp$estimate),
  HL_CI_lower = wp$conf.int[1],
  HL_CI_upper = wp$conf.int[2],
  cliff_delta = NA_real_,
  delta_CI_lower = NA_real_,
  delta_CI_upper = NA_real_,
  median_difference = NA_real_,
  median_diff_CI_lower = NA_real_,
  median_diff_CI_upper = NA_real_,
  p = wp$p.value,
  p_wilcoxon = wp$p.value
)

# 4b. Independent-group comparisons.
# Preserve the original contrast and bootstrap order.

pairs <- list(
  c(groups[1], groups[3]),
  c(groups[2], groups[3]),
  c(groups[1], groups[4]),
  c(groups[2], groups[4]),
  c(groups[3], groups[4])
)

res <- rbind(
  paired_result,
  do.call(rbind, lapply(pairs, function(k) {
    compare_groups(k[1], k[2])
  }))
)

res$p_BH <- p.adjust(res$p, method = "BH")

effect_label <- function(x) {
  as.character(cut(
    x,
    c(-Inf, .147, .33, .474, Inf),
    labels = c("negligible", "small", "medium", "large")
  ))
}

# Preserve the original reported magnitude classification.
lo <- pmin(abs(res$delta_CI_lower), abs(res$delta_CI_upper))
hi <- pmax(abs(res$delta_CI_lower), abs(res$delta_CI_upper))

res$magnitude_range <- ifelse(
  is.na(lo),
  NA_character_,
  ifelse(
    effect_label(lo) == effect_label(hi),
    effect_label(lo),
    paste(effect_label(lo), "to", effect_label(hi))
  )
)

write.csv(
  res,
  file.path(output_dir, "main_contrasts.csv"),
  row.names = FALSE
)

add_check("Paired GTEx donors", nrow(w), nrow(w) == 470)

add_check(
  "Paired Hodges-Lehmann estimate",
  round(res$HL[1], 4),
  abs(res$HL[1] - 0.2131221) < 1e-4
)

add_check(
  "Benjamini-Hochberg adjustment",
  "",
  isTRUE(all.equal(p.adjust(res$p, "BH"), res$p_BH))
)

add_check(
  "Bootstrap CIs contain point estimates",
  "",
  all(
    res$delta_CI_lower[-1] <= res$cliff_delta[-1] &
    res$cliff_delta[-1] <= res$delta_CI_upper[-1]
  )
)

add_check(
  "All six adjusted P values below 0.05",
  paste(signif(max(res$p_BH), 3), "(max)"),
  all(res$p_BH < 0.05)
)

previous_path <- file.path(
  input_root,
  "results_v2/tables/contrastes_finales_CSPG4.csv"
)

if (file.exists(previous_path)) {
  previous <- read.csv(previous_path)

  add_check(
    "Brunner-Munzel P values agree with previous analysis",
    "",
    isTRUE(all.equal(
      log10(previous$p),
      log10(res$p),
      tolerance = 1e-6
    ))
  )
}

# 5. Sensitivity analyses --------------------------------------

sens <- do.call(rbind, lapply(groups[1:3], function(g) {

  x <- v[[g]]
  y <- a65$CSPG4_log2

  data.frame(
    contrast = paste(g, "vs acral (65 primary tumors)"),
    n1 = length(x),
    n2 = length(y),
    cliff_delta = cliff_delta(x, y),
    p = brunnermunzel::brunnermunzel.test(x, y)$p.value
  )
}))

write.csv(
  sens,
  file.path(output_dir, "acral65_sensitivity.csv"),
  row.names = FALSE
)

if (requireNamespace("FSA", quietly = TRUE)) {

  dn <- suppressMessages(
    FSA::dunnTest(
      log2TPM1 ~ group,
      data = m,
      method = "bh"
    )$res
  )

  write.csv(
    dn,
    file.path(output_dir, "dunn_sensitivity.csv"),
    row.names = FALSE
  )
}

# 6. GTEx normal tissue categories -----------------------------

acral_tpm <- median(a44$CSPG4_TPM)

tissue$mayor_que_acral <- tissue$mediana_TPM > acral_tpm

tissue_names <- c(
  tejido = "tissue",
  mediana_TPM = "median_TPM",
  log2TPM1 = "log2TPM1",
  mayor_que_acral = "above_acral_median"
)

stopifnot(all(names(tissue) %in% names(tissue_names)))

names(tissue) <- unname(tissue_names[names(tissue)])

write.csv(
  tissue,
  file.path(output_dir, "normal_tissues_vs_acral.csv"),
  row.names = FALSE
)

add_check(
  "54 GTEx categories",
  nrow(tissue),
  nrow(tissue) == 54
)

add_check(
  "10 categories above the acral median",
  sum(tissue$above_acral_median),
  sum(tissue$above_acral_median) == 10
)

# 7. Optional purity and deconvolution analyses ----------------
#
# These inputs are restricted to authorized users.

supplement_path <- file.path(
  input_root, "data/raw/Supplementary_Table_1.csv"
)

if (file.exists(supplement_path)) {

  sup <- read.csv(
    supplement_path,
    check.names = FALSE
  )[, c("Sample", "Patient", "Sample_type", "Purity", "RNA_cluster")]

  sup$Purity <- suppressWarnings(as.numeric(sup$Purity))

  # Harmonize identifiers between source datasets.
  sup$sample_id <- sub("^PD", "PR", sup$Sample)

  mm <- merge(
    a[, c("sample_id", "sample_type", "en_44",
          "CSPG4_log2", "prefijo")],
    sup,
    by = "sample_id"
  )

  # Individual-level table: private output only.
  write.csv(
    mm[, c("sample_id", "Patient", "sample_type", "en_44",
           "Purity", "RNA_cluster", "CSPG4_log2")],
    file.path(private_dir, "acral_patient_purity.csv"),
    row.names = FALSE
  )

  add_check(
    "RNA samples matched to supplementary data",
    paste(nrow(mm), "of", nrow(a)),
    nrow(mm) > 0
  )

  p44 <- mm[mm$en_44, ]
  ok <- !is.na(p44$Purity)

  if (sum(ok) > 5) {

    sp <- suppressWarnings(cor.test(
      p44$CSPG4_log2[ok],
      p44$Purity[ok],
      method = "spearman"
    ))

    write.csv(
      data.frame(
        n = sum(ok),
        rho = unname(sp$estimate),
        p = sp$p.value
      ),
      file.path(output_dir, "CSPG4_vs_purity.csv"),
      row.names = FALSE
    )
  }
}

deconvolution_path <- file.path(
  input_root, "data/raw/Sample_pass_deconvolution_data.csv"
)

if (file.exists(deconvolution_path)) {

  dc <- read.csv(deconvolution_path, check.names = FALSE)

  fractions <- setdiff(names(dc), "sample")

  mm2 <- merge(
    a44[, c("sample_id", "CSPG4_log2")],
    dc,
    by.x = "sample_id",
    by.y = "sample"
  )

  ct <- do.call(rbind, lapply(fractions, function(f) {

    r <- suppressWarnings(cor.test(
      mm2$CSPG4_log2,
      mm2[[f]],
      method = "spearman"
    ))

    data.frame(
      fraction = f,
      rho = unname(r$estimate),
      p = r$p.value
    )
  }))

  ct$p_BH <- p.adjust(ct$p, "BH")

  write.csv(
    ct,
    file.path(output_dir, "CSPG4_vs_deconvolution.csv"),
    row.names = FALSE
  )

  add_check(
    "44 samples matched to deconvolution data",
    nrow(mm2),
    nrow(mm2) == 44
  )
}

# 8. Consistency checks and session record ---------------------

audit <- do.call(rbind, checks)

write.csv(
  audit,
  file.path(output_dir, "analysis_checks.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(private_dir, "sessionInfo.txt")
)

print(audit, row.names = FALSE)

cat(
  "\nChecks:", nrow(audit),
  "| Failed:", sum(!audit$ok), "\n"
)

stopifnot(all(audit$ok))

print(res, digits = 5, row.names = FALSE)

cat("\nAggregate results:", normalizePath(output_dir), "\n")
cat("Private outputs:", normalizePath(private_dir), "\n")
