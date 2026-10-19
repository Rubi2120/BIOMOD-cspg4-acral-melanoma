# ============================================================
# Funciones auxiliares (R base). Probadas con datos simulados.
# ============================================================

strip_version <- function(ids) sub("[.].*$", "", ids)   # ENSG00000141510.18 -> ENSG00000141510

write_tab <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, path, row.names = FALSE)
  message("  -> escrito: ", path)
}

# Resumen descriptivo y proporción sobre el umbral de "detectable"
summarize_expr <- function(x, detect = 1) {
  x <- x[!is.na(x)]
  data.frame(
    n = length(x),
    median = median(x), q1 = unname(quantile(x, 0.25)), q3 = unname(quantile(x, 0.75)),
    mean = mean(x), min = min(x), max = max(x),
    pct_above_detect = 100 * mean(x >= detect)
  )
}

# Tamaño de efecto de Kruskal-Wallis: epsilon^2 = H / (n - 1)  (Tomczak y Tomczak, 2014)
epsilon_sq <- function(H, n) as.numeric(H / (n - 1))

# Delta de Cliff entre dos grupos: 2*W/(n1*n2) - 1, con W de wilcox.test (x vs y)
cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  W <- suppressWarnings(wilcox.test(x, y, exact = FALSE)$statistic)
  unname(2 * W / (length(x) * length(y)) - 1)
}

# Prueba de Dunn post-hoc (Dunn, 1964) con corrección por empates y ajuste BH.
# Implementada a mano porque el paquete FSA no estaba disponible en el entorno de prueba.
# Verificada: con k = 2 grupos reproduce el p de wilcox.test(exact = FALSE, correct = FALSE).
dunn_test <- function(values, groups, p_adjust = "BH") {
  ok <- !is.na(values) & !is.na(groups)
  values <- values[ok]; groups <- droplevels(factor(groups[ok]))
  N <- length(values); r <- rank(values)
  lv <- levels(groups)
  mean_rank <- tapply(r, groups, mean); n_g <- tapply(r, groups, length)
  ties <- table(values); tie_term <- sum(ties^3 - ties) / (12 * (N - 1))
  sigma2 <- N * (N + 1) / 12 - tie_term
  out <- list()
  for (i in seq_len(length(lv) - 1)) for (j in (i + 1):length(lv)) {
    gi <- lv[i]; gj <- lv[j]
    z <- (mean_rank[gi] - mean_rank[gj]) / sqrt(sigma2 * (1 / n_g[gi] + 1 / n_g[gj]))
    out[[length(out) + 1]] <- data.frame(group1 = gi, group2 = gj, z = unname(z),
                                         p = 2 * pnorm(-abs(unname(z))))
  }
  res <- do.call(rbind, out)
  res$p_adj <- p.adjust(res$p, method = p_adjust)
  res
}

# Busca la fila de un gen por ENSG (con o sin versión) o por símbolo; devuelve índice o integer(0)
find_gene_row <- function(ids, ensg, symbol) {
  hit <- which(strip_version(ids) == ensg)
  if (length(hit) == 0) hit <- which(ids == symbol)
  hit
}
