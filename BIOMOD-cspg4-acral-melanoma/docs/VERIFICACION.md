# Lista de verificación (hacer ANTES de citar o interpretar)

Marca cada punto cuando lo hayas comprobado contra la fuente original.

## Datos y código
- [ ] Ensembl ID de CSPG4 (ENSG00000173546) confirmado en Ensembl y NCBI Gene (ID 1464).
- [ ] El enlace de Dropbox descarga y contiene los dos archivos del Nivel A.
- [ ] La matriz del Nivel A trae CSPG4 y está en escala log2(TPM+1) (el script 01 lo avisa).
- [ ] Nombres reales de columnas del CSV de metadatos (ulceración, estadio, cluster) anotados en `scripts/config.R`.
- [ ] Prefiltro de DESeq2 (≥10 counts en ≥5 muestras) confirmado en la página de protocols.io.
- [ ] Nombres exactos de los archivos de GTEx (versión) confirmados en gtexportal.org.
- [ ] Convención «detectable = TPM ≥ 1» respaldada con una cita y aceptada por el equipo.

## Citas del texto del proyecto (no verificadas)
- [ ] Ferrone et al., 2020: CSPG4 en ~80 % de los melanomas.
- [ ] Nishi et al., 2010: ~53.6 % en melanoma acral lentiginoso.
- [ ] Lino-Silva et al., 2016 (Ann Surg Oncol 23:4189-4194): melanoma acral hasta 44.1 % en México.
- [ ] Zhang et al., 2025 (LNP catiónica) y Schuler et al., 2022 (Ribociclib).

## Permisos y colaboración
- [ ] Avisar a Stef y Ferro del uso de sus archivos; definir reconocimiento (agradecimiento, cita o colaboración).
- [ ] No publicar la tabla por muestra del Nivel A hasta tener su visto bueno (está en `.gitignore`).
- [ ] Licencia del repositorio decidida (hasta entonces, no agregar archivo LICENSE).
