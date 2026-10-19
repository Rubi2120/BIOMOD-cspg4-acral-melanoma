# ============================================================
# Configuración central del proyecto
# Edita SOLO este archivo si cambian nombres de archivos o columnas.
# ============================================================

# ---- Codificación: que R lea/escriba en UTF-8 (acentos) sin importar el sistema -------
for (loc in c("C.UTF-8", "en_US.UTF-8", "es_MX.UTF-8"))
  if (!is.na(suppressWarnings(Sys.setlocale("LC_CTYPE", loc))) && nzchar(Sys.getlocale("LC_CTYPE"))) break

# ---- Gen de interés ----------------------------------------------------------
GENE_SYMBOL <- "CSPG4"
GENE_ENSG   <- "ENSG00000173546"   # PENDIENTE: confirmar en Ensembl / NCBI Gene (ID 1464)

# ---- Nivel A: melanoma acral mexicano (Basurto-Lozada et al., Nature 2026) -----
# Archivos de la carpeta de Dropbox indicada en el repositorio
# https://github.com/CGBio-Lab/Mex-acral-exomes-transcriptomes ("Input files download")
A_EXPR_FILE <- "data/raw/acralMX_expression_data_update44s.rds"
A_META_FILE <- "data/raw/acralMX_metadata_cc_update.csv"

# Nombres de columnas del CSV de metadatos.
# ¡VERIFICAR con names(meta)! Si no coinciden, el script 01 se detiene y muestra las reales.
A_COL_SAMPLE  <- "sample"
A_COLS_GROUPS <- c("ulceration", "stage")   # variables a comparar (las que no existan se omiten con aviso)
A_COL_CLUSTER <- "cluster"                  # opcional: cluster transcripcional (puede no estar en el CSV)

# ---- Nivel B: TCGA-SKCM ------------------------------------------------------
TCGA_PROJECT      <- "TCGA-SKCM"
TCGA_DIR          <- "data/raw/GDCdata"
TCGA_SAMPLE_TYPES <- c("Primary Tumor", "Metastatic")
TCGA_REFERENCE    <- "Primary Tumor"        # nivel de referencia del contraste
# Prefiltro tomado del protocolo de protocols.io (Zubair 2026): >=10 counts en >=5 muestras.
# PENDIENTE: confirmar ese renglón en la página (se ve alterado).
PREFILTER_MIN_COUNT   <- 10
PREFILTER_MIN_SAMPLES <- 5

# ---- Nivel C: GTEx (v8; descargar de gtexportal.org > Downloads) ------------------
# VERIFICAR los nombres exactos de archivo en el portal.
GTEX_TPM_GCT    <- "data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
GTEX_SAMPLES    <- "data/raw/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
GTEX_SKIN_SITES <- c("Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)")

# ---- Criterio de "detectable" -------------------------------------------------
# Convención común: TPM >= 1. En escala log2(TPM+1) equivale a >= 1.
# PENDIENTE: citar la fuente de la convención y confirmarla con el equipo.
DETECT_LOG2 <- 1

# ---- Estadística ----------------------------------------------------------------
ALPHA <- 0.05
MIN_N_PER_GROUP <- 3    # con menos muestras por grupo se avisa y no se prueba

# ---- Carpetas de salida -----------------------------------------------------------
OUT_TAB <- "results/tables"
OUT_FIG <- "results/figures"
