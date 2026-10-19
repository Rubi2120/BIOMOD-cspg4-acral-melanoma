# ============================================================
# 03 · Nivel C: CSPG4 en GTEx (piel sana y otros tejidos)
# Entrada : GTEx v8 gene TPM (.gct.gz) y SampleAttributesDS.txt (ver data/README.md)
# Salidas : results/tables/nivelC_*.csv
# Requiere: Linux/macOS (usa gzip y awk para leer solo la fila de CSPG4 sin cargar 1.5 GB en memoria).
#           En Windows, usar WSL.
# Ejecutar desde la raíz del repo:  Rscript scripts/03_nivelC_gtex.R
# ============================================================
source("scripts/config.R"); source("scripts/helpers.R")
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
for (f in c(GTEX_TPM_GCT, GTEX_SAMPLES))
  if (!file.exists(f)) stop("No existe: ", f, "\nDescárgalo de gtexportal.org (ver data/README.md) o ajusta config.R.")

# ---- 1. Leer solo encabezado (línea 3 del .gct) y la fila de CSPG4 -------------------
# Formato .gct: línea 1 = versión, línea 2 = dimensiones, línea 3 = encabezado, luego genes.
cmd <- sprintf("gzip -dc '%s' | awk -F'\\t' 'NR==3 || $1 ~ /^%s/'", GTEX_TPM_GCT, GENE_ENSG)
lines <- system(cmd, intern = TRUE)
if (length(lines) < 2) stop("No se encontró ", GENE_ENSG, " en el archivo de GTEx. Revisa el Ensembl ID.")
if (length(lines) > 2) warning("Más de una fila coincide con ", GENE_ENSG, "; se usa la primera.")
hdr  <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
vals <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
message("Gen encontrado: ", vals[1], " (", vals[2], ") | muestras en el archivo: ", length(hdr) - 2)
if (!identical(toupper(vals[2]), toupper(GENE_SYMBOL)))
  warning("El símbolo en GTEx es '", vals[2], "', no '", GENE_SYMBOL, "'. Verifica el Ensembl ID.")
tpm <- data.frame(SAMPID = hdr[-(1:2)], tpm = as.numeric(vals[-(1:2)]), stringsAsFactors = FALSE)

# ---- 2. Unir con anotación de tejido ----------------------------------------------------
att <- read.delim(GTEX_SAMPLES, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
stopifnot(all(c("SAMPID", "SMTSD") %in% names(att)))
d <- merge(tpm, att[, c("SAMPID", "SMTSD")], by = "SAMPID")
d$log2tpm1 <- log2(d$tpm + 1)
message("Muestras con TPM y anotación: ", nrow(d), " | tejidos: ", length(unique(d$SMTSD)))

# ---- 3. Piel (dos sitios, por separado) ---------------------------------------------------
skin <- d[d$SMTSD %in% GTEX_SKIN_SITES, ]
if (nrow(skin) == 0) stop("No hay muestras con los nombres de tejido de GTEX_SKIN_SITES. Tejidos disponibles que contienen 'Skin': ",
                          paste(grep("Skin", unique(d$SMTSD), value = TRUE), collapse = " | "))
write_tab(data.frame(sample = skin$SAMPID, tissue = skin$SMTSD, tpm = skin$tpm, log2tpm1 = skin$log2tpm1),
          file.path(OUT_TAB, "nivelC_gtex_piel_por_muestra.csv"))
print(do.call(rbind, lapply(split(skin, skin$SMTSD),
                            function(x) cbind(tejido = x$SMTSD[1], summarize_expr(x$log2tpm1, DETECT_LOG2)))))

# ---- 4. Todos los tejidos (objetivo 4: especificidad / riesgo fuera del blanco) ---------------
all_tis <- do.call(rbind, lapply(split(d, d$SMTSD), function(x)
  data.frame(tejido = x$SMTSD[1], n = nrow(x), mediana_tpm = median(x$tpm), q1_tpm = unname(quantile(x$tpm, .25)),
             q3_tpm = unname(quantile(x$tpm, .75)))))
all_tis <- all_tis[order(-all_tis$mediana_tpm), ]
write_tab(all_tis, file.path(OUT_TAB, "nivelC_gtex_todos_tejidos.csv"))
message("Top 5 tejidos por mediana de TPM:"); print(head(all_tis, 5), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_TAB, "nivelC_sessionInfo.txt"))
