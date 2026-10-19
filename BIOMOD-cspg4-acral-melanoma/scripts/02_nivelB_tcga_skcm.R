# ============================================================
# 02 · Nivel B: CSPG4 en TCGA-SKCM (primario vs metástasis) con DESeq2
# Flujo adaptado de: Zubair H. (2026) protocols.io, doi:10.17504/protocols.io.kqdg3mn6el25/v1
# (validado por su autor en TCGA-LUAD; los valores de referencia que reporta NO aplican a SKCM).
# Referencias primarias: Love et al. 2014 (DESeq2), Colaprico et al. 2016 (TCGAbiolinks).
#
# ESTADO: NO PROBADO en el entorno de desarrollo (requiere descargar ~GB del GDC).
#         Revisa cada mensaje y detén el script si algo no coincide.
# Requiere: R >= 4.4 con Bioconductor (TCGAbiolinks, SummarizedExperiment, DESeq2, apeglm), ~20 GB libres.
# Ejecutar desde la raíz del repo:  Rscript scripts/02_nivelB_tcga_skcm.R
# ============================================================
source("scripts/config.R"); source("scripts/helpers.R")
suppressPackageStartupMessages({ library(TCGAbiolinks); library(SummarizedExperiment); library(DESeq2) })
set.seed(42)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE); dir.create(TCGA_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Descargar counts crudos (flujo STAR - Counts del GDC) ---------------------------------
query <- GDCquery(project = TCGA_PROJECT, data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(query, directory = TCGA_DIR)   # si falla a la mitad: borrar la carpeta parcial y repetir (no reanudar)
se <- GDCprepare(query, directory = TCGA_DIR)

message("Assays disponibles: ", paste(assayNames(se), collapse = ", "))
stopifnot(all(c("unstranded", "tpm_unstrand") %in% assayNames(se)))   # si falla, ver assayNames(se)
message("Tipos de muestra en TCGA-SKCM:"); print(table(colData(se)$sample_type))

# ---- 2. Quedarse con primario y metástasis ------------------------------------------------------
keep_s <- colData(se)$sample_type %in% TCGA_SAMPLE_TYPES
se <- se[, keep_s]
message("Muestras usadas: ", ncol(se)); print(table(colData(se)$sample_type))
# Nota: algunos pacientes aportan más de una muestra; se asume independencia (limitación a declarar).

# ---- 3. DESeq2 con COUNTS crudos (nunca TPM/FPKM) -----------------------------------------------
counts_mat <- assay(se, "unstranded")
col_data   <- colData(se)
stopifnot(identical(colnames(counts_mat), rownames(col_data)))
dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~ sample_type)
keep <- rowSums(counts(dds) >= PREFILTER_MIN_COUNT) >= PREFILTER_MIN_SAMPLES
dds  <- dds[keep, ]
dds$sample_type <- relevel(factor(dds$sample_type), ref = TCGA_REFERENCE)   # el primer nivel es la referencia
message("Niveles (el primero es la referencia): ", paste(levels(dds$sample_type), collapse = " | "))
dds <- DESeq(dds)
message("Coeficientes: ", paste(resultsNames(dds), collapse = ", "))
res <- results(dds, name = resultsNames(dds)[2], alpha = ALPHA)
res_shr <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")

# ---- 4. Fila de CSPG4 -------------------------------------------------------------------------------------
ids <- strip_version(rownames(dds))
i <- which(ids == GENE_ENSG)
if (length(i) == 0) {
  message("CSPG4 no pasó el prefiltro o el ID no coincide; revisar. Se continúa solo con los valores por muestra.")
  de_row <- NULL
} else {
  de_row <- data.frame(gene = GENE_SYMBOL, ensg = GENE_ENSG,
                       comparacion = paste(resultsNames(dds)[2], "(referencia:", TCGA_REFERENCE, ")"),
                       baseMean = res$baseMean[i], log2FC = res$log2FoldChange[i],
                       log2FC_apeglm = res_shr$log2FoldChange[i], padj = res$padj[i],
                       n_primario = sum(dds$sample_type == "Primary Tumor"),
                       n_metastasis = sum(dds$sample_type == "Metastatic"))
  print(de_row); write_tab(de_row, file.path(OUT_TAB, "nivelB_deseq2_cspg4.csv"))
}
# Umbrales del protocolo: padj < 0.05 y |log2FC| > 1 (criterio estándar de DEG)

# ---- 5. TPM por muestra de CSPG4 (solo para visualizar e integrar con los otros niveles) ------------
rd_ids <- strip_version(rownames(se)); j <- which(rd_ids == GENE_ENSG)
stopifnot(length(j) >= 1)
tpm <- assay(se, "tpm_unstrand")[j[1], ]
out <- data.frame(sample = colnames(se), sample_type = colData(se)$sample_type,
                  tpm = as.numeric(tpm), log2tpm1 = log2(as.numeric(tpm) + 1))
write_tab(out, file.path(OUT_TAB, "nivelB_tcga_por_muestra.csv"))
print(do.call(rbind, lapply(split(out, out$sample_type), function(x)
  cbind(tipo = x$sample_type[1], summarize_expr(x$log2tpm1, DETECT_LOG2)))))
writeLines(c(paste("GDC / TCGAbiolinks:", as.character(packageVersion("TCGAbiolinks"))),
             capture.output(sessionInfo())), file.path(OUT_TAB, "nivelB_sessionInfo.txt"))
