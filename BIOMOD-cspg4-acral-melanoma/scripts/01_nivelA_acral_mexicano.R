# ============================================================
# 01 · Nivel A: CSPG4 en melanoma acral mexicano (44 primarios)
# Fuente: Basurto-Lozada et al., Nature 2026 (matriz y metadatos del repositorio de los autores)
# Entrada : data/raw/acralMX_expression_data_update44s.rds y data/raw/acralMX_metadata_cc_update.csv
# Salidas : results/tables/nivelA_*.csv, results/figures/nivelA_*.png
# Ejecutar desde la raíz del repo:  Rscript scripts/01_nivelA_acral_mexicano.R
# ============================================================
source("scripts/config.R"); source("scripts/helpers.R")
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE); dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Cargar y revisar la matriz -----------------------------------------------
for (f in c(A_EXPR_FILE, A_META_FILE))
  if (!file.exists(f)) stop("No existe: ", f, "\nDescárgalo de la carpeta de Dropbox (ver data/README.md).")

mt <- readRDS(A_EXPR_FILE)
mt <- as.matrix(mt)
if (!is.numeric(mt)) stop("La matriz no es numérica; revisa el objeto con str(readRDS(...)).")
message("Matriz: ", nrow(mt), " filas x ", ncol(mt), " columnas")
message("Primeros nombres de fila: ", paste(head(rownames(mt), 3), collapse = ", "))

# Orientación: se espera genes en filas y muestras en columnas
if (length(find_gene_row(rownames(mt), GENE_ENSG, GENE_SYMBOL)) == 0 &&
    length(find_gene_row(colnames(mt), GENE_ENSG, GENE_SYMBOL)) > 0) {
  message("Aviso: la matriz parece estar transpuesta; se transpone.")
  mt <- t(mt)
}

# Escala: log2(TPM+1) no puede superar ~20 (TPM <= 1e6). Si pasa de 21, es TPM lineal.
mx <- max(mt, na.rm = TRUE)
message("Valor máximo de la matriz: ", round(mx, 2))
if (mx > 21) {
  message("Aviso: máximo > 21, parece TPM lineal. Se transforma a log2(TPM+1).")
  mt <- log2(mt + 1)
} else {
  message("Escala compatible con log2(TPM+1) (coincide con lo descrito en los Methods del artículo).")
}

# ---- 2. Extraer CSPG4 -------------------------------------------------------------------
idx <- find_gene_row(rownames(mt), GENE_ENSG, GENE_SYMBOL)
if (length(idx) == 0)
  stop(GENE_SYMBOL, " NO está en la matriz (posible filtrado de genes). ",
       "Pedir la matriz completa a las autoras/autores o tramitar acceso a EGA (EGAD00001015756).")
if (length(idx) > 1) { message("Aviso: ", length(idx), " filas coinciden; se usa la primera."); idx <- idx[1] }
cspg4 <- data.frame(sample = colnames(mt), log2tpm1 = as.numeric(mt[idx, ]), stringsAsFactors = FALSE)

# ---- 3. Unir con metadatos ----------------------------------------------------------------
meta <- read.csv(A_META_FILE, stringsAsFactors = FALSE, check.names = FALSE)
message("Columnas del CSV de metadatos: ", paste(names(meta), collapse = " | "))
if (!A_COL_SAMPLE %in% names(meta))
  stop("La columna '", A_COL_SAMPLE, "' no existe en los metadatos. Ajusta A_COL_SAMPLE en scripts/config.R ",
       "usando una de las columnas listadas arriba.")
d <- merge(cspg4, meta, by.x = "sample", by.y = A_COL_SAMPLE, all.x = TRUE)
n_sin_meta <- sum(!cspg4$sample %in% meta[[A_COL_SAMPLE]])
message("Muestras en la matriz: ", nrow(cspg4), " | con metadatos: ", nrow(cspg4) - n_sin_meta)
if (n_sin_meta > 0) warning(n_sin_meta, " muestras sin metadatos: revisa que los identificadores coincidan.")
write_tab(d, file.path(OUT_TAB, "nivelA_cspg4_por_muestra.csv"))

# ---- 4. Descriptivo general (distribución, no solo media) ----------------------------------
desc <- cbind(grupo = "Todos (acral MX)", summarize_expr(d$log2tpm1, DETECT_LOG2))
print(desc)

p_dist <- ggplot(d, aes(log2tpm1)) +
  geom_histogram(bins = 15, fill = "grey70", colour = "white") +
  geom_vline(xintercept = DETECT_LOG2, linetype = "dashed") +
  labs(x = paste0(GENE_SYMBOL, " (log2(TPM+1))"), y = "Número de tumores",
       title = "Distribución de CSPG4 en melanoma acral mexicano",
       subtitle = "Línea punteada: umbral de «detectable» (TPM ≥ 1; convención por confirmar)") +
  theme_bw()
ggsave(file.path(OUT_FIG, "nivelA_distribucion.png"), p_dist, width = 6, height = 4, dpi = 300)

# ---- 5. Comparaciones pre-especificadas ---------------------------------------------------------
vars <- A_COLS_GROUPS
if (A_COL_CLUSTER %in% names(d)) vars <- c(vars, A_COL_CLUSTER)
faltan <- setdiff(vars, names(d))
if (length(faltan)) message("Variables no encontradas (se omiten): ", paste(faltan, collapse = ", "))
vars <- intersect(vars, names(d))

res_list <- list(); desc_list <- list(desc)
for (v in vars) {
  g <- factor(d[[v]]); keep <- !is.na(g); dd <- d[keep, ]; dd$grupo <- droplevels(g[keep])
  tab <- table(dd$grupo)
  message("\n== ", v, " ==  n por grupo: ", paste(names(tab), tab, sep = "=", collapse = ", "))
  for (lv in names(tab))
    desc_list[[length(desc_list) + 1]] <- cbind(grupo = paste0(v, ": ", lv),
                                                summarize_expr(dd$log2tpm1[dd$grupo == lv], DETECT_LOG2))
  if (nlevels(dd$grupo) < 2 || any(tab < MIN_N_PER_GROUP)) {
    message("  Se omite la prueba: menos de ", MIN_N_PER_GROUP, " muestras en algún grupo.")
    next
  }
  if (nlevels(dd$grupo) == 2) {
    lv <- levels(dd$grupo); x <- dd$log2tpm1[dd$grupo == lv[1]]; y <- dd$log2tpm1[dd$grupo == lv[2]]
    wt <- wilcox.test(x, y, exact = FALSE)
    res_list[[v]] <- data.frame(variable = v, prueba = "Wilcoxon-Mann-Whitney (2 colas)",
                                comparacion = paste(lv, collapse = " vs "), estadistico = unname(wt$statistic),
                                p = wt$p.value, efecto = cliffs_delta(x, y), tipo_efecto = "delta de Cliff")
  } else {
    kw <- kruskal.test(log2tpm1 ~ grupo, data = dd)
    res_list[[v]] <- data.frame(variable = v, prueba = "Kruskal-Wallis",
                                comparacion = paste(levels(dd$grupo), collapse = " / "),
                                estadistico = unname(kw$statistic), p = kw$p.value,
                                efecto = epsilon_sq(unname(kw$statistic), nrow(dd)), tipo_efecto = "epsilon cuadrado")
  }
  pl <- ggplot(dd, aes(grupo, log2tpm1)) +
    geom_boxplot(outlier.shape = NA, fill = NA) + geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
    geom_hline(yintercept = DETECT_LOG2, linetype = "dashed") +
    labs(x = v, y = paste0(GENE_SYMBOL, " (log2(TPM+1))"),
         title = paste("CSPG4 por", v, "(acral mexicano)")) + theme_bw()
  ggsave(file.path(OUT_FIG, paste0("nivelA_", v, ".png")), pl, width = 5, height = 4, dpi = 300)
}

write_tab(do.call(rbind, desc_list), file.path(OUT_TAB, "nivelA_descriptivo.csv"))
if (length(res_list)) {
  res <- do.call(rbind, res_list); rownames(res) <- NULL
  res$p_ajustado_BH <- p.adjust(res$p, method = "BH")   # ajuste por las comparaciones pre-especificadas
  print(res); write_tab(res, file.path(OUT_TAB, "nivelA_pruebas.csv"))
}
message("\nListo. Revisa results/tables y results/figures. Anota versiones con sessionInfo() en docs/.")
writeLines(capture.output(sessionInfo()), file.path(OUT_TAB, "nivelA_sessionInfo.txt"))
