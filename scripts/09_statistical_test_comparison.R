# BIOMOD 2026 - Pandora's Box
# Comparison of statistical procedures for CSPG4 expression.
#
# Compare Dunn, Wilcoxon rank-sum and Brunner-Munzel tests
# using the same samples as the main expression analysis.
#
# Dunn uses pooled ranks across four groups and BH adjustment
# across six comparisons. The other tests use pairwise data.
#
# The omnibus Kruskal-Wallis and Dunn analyses are retained
# as exploratory sensitivity analyses. They do not account
# for pairing between the two GTEx skin groups.

local({

  suppressPackageStartupMessages(library(brunnermunzel))

  input_root <- getOption(
    "biomod.tests.input_root", "."
  )

  output_dir <- getOption(
    "biomod.tests.output_dir", "tables"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

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

  non_sun <- grepl("Not Sun Exposed", p$SMTSD)

  v <- c(
    p$CSPG4_log2TPM1[non_sun],
    p$CSPG4_log2TPM1[!non_sun],
    b$CSPG4_log2TPM1[b$sample_type == "Primary Tumor"],
    a$CSPG4_log2[a$en_44]
  )

  grp <- factor(
    rep(
      c("G1a", "G1b", "G2", "G3"),
      c(
        sum(non_sun),
        sum(!non_sun),
        sum(b$sample_type == "Primary Tumor"),
        sum(a$en_44)
      )
    ),
    levels = c("G1a", "G1b", "G2", "G3")
  )

  print(table(grp))

  stopifnot(identical(
    as.numeric(table(grp)),
    c(604, 701, 103, 44)
  ))

  # Dunn's test with correction for tied ranks.

  dunn_bh <- function(v, g) {

    N <- length(v)
    r <- rank(v)
    tt <- table(v)

    s2 <- N * (N + 1) / 12 -
      sum(tt^3 - tt) / (12 * (N - 1))

    mean_ranks <- tapply(r, g, mean)

    group_sizes <- table(g)

    out <- do.call(rbind, lapply(
      combn(levels(g), 2, simplify = FALSE),
      function(q) {

        z <- (
          mean_ranks[[q[1]]] - mean_ranks[[q[2]]]
        ) / sqrt(
          s2 * (
            1 / group_sizes[[q[1]]] +
            1 / group_sizes[[q[2]]]
          )
        )

        data.frame(
          a = q[1],
          b = q[2],
          p_raw = 2 * pnorm(-abs(z))
        )
      }
    ))

    out$p_BH <- p.adjust(out$p_raw, "BH")

    out
  }

  # Omnibus comparison of the four groups.

  kw <- kruskal.test(v ~ grp)

  cat(sprintf(
    "\nKruskal-Wallis: chi-square = %.1f, p = %.2e\n",
    kw$statistic,
    kw$p.value
  ))

  dn <- dunn_bh(v, grp)

  # The five independent-group contrasts.
  # Dunn's adjustment retains the original six-test family.

  pairs <- list(
    c("G1a", "G2"),
    c("G1b", "G2"),
    c("G1a", "G3"),
    c("G1b", "G3"),
    c("G2", "G3")
  )

  results <- do.call(rbind, lapply(pairs, function(q) {

    x <- v[grp == q[1]]
    y <- v[grp == q[2]]

    d <- dn[
      dn$a == q[1] & dn$b == q[2],
    ]

    data.frame(
      contrast = paste(q, collapse = " vs "),
      p_Dunn_raw = d$p_raw,
      p_Dunn_BH = d$p_BH,
      p_Wilcoxon_raw = wilcox.test(
        x, y, exact = FALSE
      )$p.value,
      p_BrunnerMunzel_raw =
        brunnermunzel.test(x, y)$p.value
    )
  }))

  stopifnot(
    nrow(dn) == 6,
    nrow(results) == 5
  )

  print(results, digits = 6, row.names = FALSE)

  write.csv(
    results,
    file.path(
      output_dir,
      "statistical_test_comparison.csv"
    ),
    row.names = FALSE
  )

  # Interpretation:
  # The procedures differ in their rank calculations,
  # null hypotheses and multiple-testing adjustments.
  # Their P values should not be interpreted as directly
  # equivalent measures of evidence.
  #
  # Dunn's G2 versus G3 comparison is not significant at
  # the 0.05 threshold, whereas the pairwise tests yield
  # smaller P values.

  cat("\nG2 versus G3:\n")

  print(
    results[results$contrast == "G2 vs G3", ],
    digits = 6,
    row.names = FALSE
  )

})
