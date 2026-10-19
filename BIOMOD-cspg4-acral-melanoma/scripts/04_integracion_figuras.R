# ============================================================
# 04 · Integración EXPLORATORIA de los tres niveles + figura principal
# Entradas: tablas por muestra de los scripts 01, 02 y 03 (results/tables/).
# AVISO: la comparación entre fuentes es exploratoria. Cada grupo viene de una sola fuente, por lo que
# fuente y condición están confundidas; NO se aplica ComBat/removeBatchEffect (ver protocolo, sec. 5.4).
# Ejecutar desde la raíz del repo:  Rscript scripts/04_integracion_figuras.R
# ============================================================
source("scripts/config.R"); source("scripts/helpers.R")
suppressPackageStartupMessages(library(ggplot2))

fA <- file.path(OUT_TAB, "nivelA_cspg4_por_muestra.csv")
fB <- file.path(OUT_TAB, "nivelB_tcga_por_muestra.csv")
fC <- file.path(OUT_TAB, "nivelC_gtex_piel_por_muestra.csv")
falt <- c(fA, fB, fC)[!file.exists(c(fA, fB, fC))]
if (length(falt)) stop("Faltan tablas de entrada (corre antes los scripts 01-03): ", paste(falt, collapse = ", "))

a <- read.csv(fA); b <- read.csv(fB); cc <- read.csv(fC)
d <- rbind(
  data.frame(fuente = "GTEx", grupo = ifelse(grepl("Suprapubic", cc$tissue), "Piel GTEx (no expuesta)", "Piel GTEx (expuesta)"),
             sample = cc$sample, log2tpm1 = cc$log2tpm1),
  data.frame(fuente = "TCGA-SKCM", grupo = ifelse(b$sample_type == "Primary Tumor", "SKCM primario (TCGA)", "SKCM metástasis (TCGA)"),
             sample = b$sample, log2tpm1 = b$log2tpm1),
  data.frame(fuente = "Acral MX", grupo = "Acral primario (MX)", sample = a$sample, log2tpm1 = a$log2tpm1)
)
orden <- c("Piel GTEx (no expuesta)", "Piel GTEx (expuesta)", "SKCM primario (TCGA)", "SKCM metástasis (TCGA)", "Acral primario (MX)")
d$grupo <- factor(d$grupo, levels = intersect(orden, unique(d$grupo)))
write_tab(d, file.path(OUT_TAB, "integrado_cspg4_log2tpm1.csv"))

# ---- Descriptivo por grupo ----
desc <- do.call(rbind, lapply(levels(d$grupo), function(g) cbind(grupo = g, summarize_expr(d$log2tpm1[d$grupo == g], DETECT_LOG2))))
print(desc); write_tab(desc, file.path(OUT_TAB, "integrado_descriptivo.csv"))

# ---- Kruskal-Wallis + Dunn (BH) + tamaños de efecto ----
kw <- kruskal.test(log2tpm1 ~ grupo, data = d)
kw_tab <- data.frame(prueba = "Kruskal-Wallis", H = unname(kw$statistic), gl = unname(kw$parameter), p = kw$p.value,
                     epsilon2 = epsilon_sq(unname(kw$statistic), nrow(d)), n_total = nrow(d), nota = "EXPLORATORIO (inter-estudio, sin corrección de batch)")
print(kw_tab); write_tab(kw_tab, file.path(OUT_TAB, "integrado_kruskal.csv"))
dn <- dunn_test(d$log2tpm1, d$grupo, "BH")
dn$delta_cliff <- mapply(function(g1, g2) cliffs_delta(d$log2tpm1[d$grupo == g1], d$log2tpm1[d$grupo == g2]), dn$group1, dn$group2)
print(dn); write_tab(dn, file.path(OUT_TAB, "integrado_dunn_BH.csv"))

# ---- Figura principal: violín + boxplot + puntos ----
p <- ggplot(d, aes(grupo, log2tpm1, fill = fuente)) +
  geom_violin(alpha = 0.35, colour = NA, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.08, height = 0, size = 0.8, alpha = 0.4) +
  geom_hline(yintercept = DETECT_LOG2, linetype = "dashed") +
  labs(x = NULL, y = paste0(GENE_SYMBOL, " (log2(TPM+1))"), fill = "Fuente",
       title = "Expresión de CSPG4 en piel sana y melanoma",
       subtitle = "Comparación exploratoria entre estudios (sin corrección de batch). Línea punteada: TPM ≥ 1") +
  scale_x_discrete(labels = function(x) sub(" [(]", "\n(", x)) + theme_bw()
ggsave(file.path(OUT_FIG, "fig1_cspg4_tres_fuentes.png"), p, width = 8, height = 5, dpi = 300)
message("Listo: results/figures/fig1_cspg4_tres_fuentes.png")
