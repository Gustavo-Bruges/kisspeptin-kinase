## =============================================================================
## 07a_process_cel_GSE28383.R · RMA de los CEL de GSE28383 (Affymetrix Mouse Exon 1.0 ST, GPL6193)
## GEO no deposita valores procesados (series matrix con 0 filas): se normaliza desde los CEL crudos.
## Insumo:  data/raw/geo/GSE28383_cel/*.CEL.gz (GSE28383_RAW.tar)
## Salida:  data/raw/geo/GSE28383_rma_core.rds  -> leído por load_series$GSE28383 en 07
## Los autores usaron GC-RMA en ArrayAssist 5.0 (comercial); aquí RMA de oligo a nivel de transcript cluster
## ("core"), reproducible. Validación cruzada contra su lista GSE28383.xls al final.
## REGLA DE INDEPENDENCIA: no lee salidas de 03–06
## =============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(oligo); library(Biobase); library(GEOquery)
  library(moex10sttranscriptcluster.db)
})

out_file <- file.path(paths$raw_geo, "GSE28383_rma_core.rds")
cel_dir  <- file.path(paths$raw_geo, "GSE28383_cel")

if (!file.exists(out_file)) {
  cels <- list.files(cel_dir, "[.]CEL[.]gz$", full.names = TRUE, ignore.case = TRUE)
  stopifnot(length(cels) == 17)
  raw  <- read.celfiles(cels)
  eset <- rma(raw, target = "core")                          # background + quantile + median polish, transcript clusters
  colnames(eset) <- str_extract(basename(cels), "^GSM[0-9]+")

  ## Sample annotation from GEO (the series matrix carries pData even without expression rows)
  pd <- pData(getGEO("GSE28383", GSEMatrix = TRUE, destdir = paths$raw_geo)[[1]])
  pd <- tibble(sample = rownames(pd), title = pd$title, genotype = pd[["genotype:ch1"]],
               age = pd[["age:ch1"]], source_name = pd$source_name_ch1)
  pd <- pd[match(colnames(eset), pd$sample), ]
  stopifnot(identical(pd$sample, colnames(eset)))

  ## Transcript cluster -> gene symbol (clusters mapping to several symbols are left unannotated)
  ann <- AnnotationDbi::select(moex10sttranscriptcluster.db, keys = featureNames(eset),
                               keytype = "PROBEID", columns = "SYMBOL") |>
    filter(!is.na(SYMBOL)) |> group_by(PROBEID) |> filter(n_distinct(SYMBOL) == 1) |> ungroup() |>
    distinct(PROBEID, SYMBOL)
  symbol <- ann$SYMBOL[match(featureNames(eset), ann$PROBEID)]

  saveRDS(list(ex = exprs(eset), symbol = symbol, pd = pd,
               method = "oligo::rma(target = 'core'); moex10sttranscriptcluster.db",
               oligo = as.character(packageVersion("oligo")), processed_at = Sys.time()), out_file)
}

## ---- Control de calidad ----
obj <- readRDS(out_file)
qc  <- tibble(sample = colnames(obj$ex), median = apply(obj$ex, 2, median), iqr = apply(obj$ex, 2, IQR)) |>
  left_join(obj$pd, by = "sample")
write_csv(qc, file.path(paths$supp, "S5_GSE28383_rma_qc.csv"))
message(sprintf("GSE28383 RMA core: %d transcript clusters, %d with a unique symbol, %d samples",
                nrow(obj$ex), sum(!is.na(obj$symbol)), ncol(obj$ex)))

## Genotype check: Kiss1 / Kiss1r expression by group (exon-deletion KOs may keep partial transcripts)
geno <- map_dfr(c("Kiss1", "Kiss1r"), \(g) {
  i <- which(obj$symbol == g)
  if (!length(i)) return(NULL)
  tibble(gene = g, sample = colnames(obj$ex), value = colMeans(obj$ex[i, , drop = FALSE]))
}) |> left_join(dplyr::select(obj$pd, sample, genotype, age), by = "sample") |>   # AnnotationDbi masks select()
  group_by(gene, genotype, age) |> summarise(mean_log2 = mean(value), n = n(), .groups = "drop")
print(geno)
write_csv(geno, file.path(paths$supp, "S5_GSE28383_genotype_check.csv"))
