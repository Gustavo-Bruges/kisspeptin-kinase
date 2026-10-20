## =============================================================================
## 04_flags_context.R · Marcas S (quinasas), E (efectores), contexto oncológico
## S y E son MARCAS sobre el corpus C, nunca filtros de búsqueda
## =============================================================================
source(here::here("R", "00_setup.R"))

M_main    <- load_versioned(paste0("M_", params$main_base))   # same base as in 03 (docs/decisions.md)
records   <- load_versioned("records_doi_pmid")
comm_main <- load_versioned("communities_main")

pm_main <- if ("PM" %in% names(M_main)) na_if(str_trim(as.character(M_main$PM)), "") else NA_character_
corpus <- M_main |>
  transmute(key = record_key(DI, UT), doi = DI, pm_native = pm_main,
            year = PY, title = TI, abstract = AB, doc_type = DT) |>
  left_join(select(records, key, pmid_resolved = pmid), by = "key") |>
  mutate(pmid = coalesce(pm_native, pmid_resolved)) |>          # native PM first; DOI lookup from 02 as fallback
  select(-pm_native, -pmid_resolved) |>
  left_join(select(comm_main, doi, community), by = "doi") |>
  mutate(
    community   = coalesce(community, "unassigned"),
    text_ta     = str_c(coalesce(title, ""), " ", coalesce(abstract, "")),
    no_abstract = is.na(abstract) | abstract == "" |
                  str_detect(abstract, regex("^\\[no abstract", ignore_case = TRUE)),   # Scopus placeholder "[No abstract available]"
    is_review   = str_detect(str_to_lower(doc_type), "review")
  )

## ---- Bloque K v1.1 (quinasas y fosforilación) ----
kinase_regex <- regex(paste0(
  "phosphorylat|kinase|phosphoproteom|kinom|",
  "\\bMAPK|\\bERK[0-9/]*\\b|\\bp38\\b|\\bJNK|\\bMEK[0-9]?\\b|mitogen-activated protein kinase|",
  "\\bPI3K|phosphatidylinositol 3-kinase|\\bAKT[0-9]?\\b|protein kinase B|\\bmTOR|",
  "\\bPKC|protein kinase C|\\bCaMK|calmodulin-dependent kinase|\\bPKA\\b|protein kinase A|",
  "\\bGRK[0-9]?\\b|G protein-coupled receptor kinase|arrestin|",
  "\\bGSK-?3|glycogen synthase kinase|\\bc-Src\\b|Src kinase|Src family|\\bFAK\\b|focal adhesion kinase|\\bp?90?RSK|",
  "\\bEIF2AK2\\b|protein kinase R|double-stranded RNA-dependent protein kinase|\\bPKR\\b|",
  "\\bAXL\\b|\\bRTKs?\\b|\\bEGFR\\b|\\bErbB[0-9]?\\b|epidermal growth factor receptor|",
  "\\bTGFBR[12]?\\b|\\bALK5\\b|TGF-(beta|\u03b2) receptor|",
  "\\bROCK[12]\\b|Rho-associated kinase|Rho kinase"            # no bare ROCK: avoids "rockfish" in teleost papers
), ignore_case = TRUE)

## ---- Bloque E (efectores no quinasa de la cascada) ----
effector_regex <- regex(paste0(
  "phospholipase C|\\bPLC[\u03b2\u03b3]?\\b|\\bG(alpha)?q(/11)?\\b|G\u03b1q|\\bRhoA\\b|",
  "\\bMMP-?[0-9]+|matrix metalloproteinase|\\bCXCR4\\b|\\bSDF-1\\b|\\bCXCL12\\b|",
  "\\bBCRP\\b|\\bABCG2\\b|invadopodi"
), ignore_case = TRUE)

corpus <- corpus |>
  mutate(signaling = str_detect(text_ta, kinase_regex),
         effector  = str_detect(text_ta, effector_regex),
         kblock_v  = params$query_version)

## ---- Contexto oncológico ----
ctx_dict <- c(
  breast      = "breast (cancer|carcinoma|tumou?r)",
  tnbc        = "triple[- ]negative|\\bTNBC\\b",
  er_neg      = "\\bER[\u03b1-]?[- ]?negative|estrogen receptor[- ](alpha[- ])?negative|\\bER\\(-\\)",
  er_pos      = "\\bER[\u03b1-]?[- ]?positive|estrogen receptor[- ](alpha[- ])?positive|\\bER\\(\\+\\)",
  colorectal  = "colorectal|colon (cancer|carcinoma)|rectal cancer|\\bCRC\\b",
  endometrial = "endometrial (cancer|carcinoma|adenocarcinoma)",   # excludes endometriosis by construction
  melanoma    = "melanoma",
  prostate    = "prostat\\w* (cancer|carcinoma|adenocarcinoma)",
  lung        = "lung (cancer|carcinoma|adenocarcinoma)|\\bNSCLC\\b",
  glioma      = "gliom|glioblastoma",
  cervical    = "cervical (cancer|carcinoma)",
  gastric     = "gastric (cancer|carcinoma)",
  ovarian     = "ovarian (cancer|carcinoma)",
  liver       = "hepatocellular|\\bHCC\\b",
  pancreatic  = "pancreatic (cancer|adenocarcinoma|ductal)"
)
ctx_tags <- map_dfc(ctx_dict, \(rx) str_detect(corpus$text_ta, regex(rx, ignore_case = TRUE)))
corpus <- bind_cols(corpus, ctx_tags) |>
  mutate(n_ctx = rowSums(across(all_of(names(ctx_dict)))))      # n_ctx > 1: multi-cancer paper, enters all its strata

## ---- Resumen S/C por año (excluye registros sin resumen en numerador y denominador) ----
signal_trend <- corpus |>
  filter(!no_abstract) |>
  count(year, community, signaling) |>
  group_by(year, community) |>
  mutate(prop = n / sum(n)) |>
  ungroup()
write_csv(signal_trend, file.path(paths$tables, "signaling_share_by_year.csv"))

## ---- Plantilla de anotación manual de dirección ((oncología ∪ molecular) ∩ S) ----
# Molecular = receptor discovery/pharmacology community (cluster 19, decisions.md 2026-10-18), included by GB
set.seed(params$seed)
template <- corpus |>
  filter(community %in% c("oncology", "molecular"), signaling, !no_abstract) |>
  mutate(direction = NA_character_,   # suppressive / promotive / mixed / not_stated
         model     = NA_character_,   # cell_line / xenograft / patient_tissue / review / other
         annotator = NA_character_,
         notes     = NA_character_) |>
  select(pmid, doi, year, community, title, abstract, all_of(names(ctx_dict)), direction, model, annotator, notes) |>
  slice_sample(prop = 1)              # random order avoids annotation drift by year or journal
write_csv(template, file.path(paths$annotation, "direction_annotation_template.csv"))

## ---- Una plantilla .xlsx por anotador (Excel no altera separador, fechas ni DOI) ----
# Dropdown lists constrain direction and model to the coded values; an existing file is never overwritten,
# so re-running 04 cannot erase an annotation in progress
direction_levels <- c("suppressive", "promotive", "mixed", "not_stated")
model_levels     <- c("cell_line", "xenograft", "patient_tissue", "review", "other")

write_annotator_xlsx <- function(df, annotator_id) {
  f <- file.path(paths$annotation, sprintf("direction_annotator%d.xlsx", annotator_id))
  if (file.exists(f)) { message("kept existing ", basename(f)); return(invisible(f)) }
  df$annotator <- sprintf("annotator_%d", annotator_id)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "anotacion")
  openxlsx::addWorksheet(wb, "instrucciones")
  openxlsx::addWorksheet(wb, "listas", visible = FALSE)
  openxlsx::writeData(wb, "listas", data.frame(direction = c(direction_levels, NA), model = model_levels))
  openxlsx::writeData(wb, "anotacion", df, withFilter = TRUE,
                      headerStyle = openxlsx::createStyle(textDecoration = "bold", fgFill = "#DDEBF7"))
  rows <- seq_len(nrow(df)) + 1
  col  <- function(x) match(x, names(df))
  openxlsx::dataValidation(wb, "anotacion", cols = col("direction"), rows = rows, type = "list",
                           value = sprintf("'listas'!$A$2:$A$%d", length(direction_levels) + 1))
  openxlsx::dataValidation(wb, "anotacion", cols = col("model"), rows = rows, type = "list",
                           value = sprintf("'listas'!$B$2:$B$%d", length(model_levels) + 1))
  openxlsx::addStyle(wb, "anotacion", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                     rows = rows, cols = c(col("title"), col("abstract")), gridExpand = TRUE)
  openxlsx::setColWidths(wb, "anotacion", cols = c(col("title"), col("abstract")), widths = c(45, 90))
  openxlsx::freezePane(wb, "anotacion", firstRow = TRUE)
  openxlsx::writeData(wb, "instrucciones", data.frame(campo = c(
    "Anotador", "Tarea", "direction", "", "", "", "model", "notes", "Regla"), descripcion = c(
    sprintf("annotator_%d (anotar sin consultar al otro anotador)", annotator_id),
    "Leer título y resumen; decidir el efecto de KISS1/KISS1R sobre el fenotipo tumoral descrito",
    "suppressive: KISS1/KISS1R inhibe invasión, migración, metástasis o crecimiento",
    "promotive: KISS1/KISS1R favorece invasión, migración, metástasis o crecimiento",
    "mixed: el artículo informa ambos efectos (p. ej., según subtipo o estado ER)",
    "not_stated: el resumen no permite decidir (incluye artículos no oncológicos)",
    "Sistema principal del estudio: cell_line, xenograft, patient_tissue, review, other",
    "Texto libre: dudas, frase del resumen que justifica la decisión",
    "No modificar pmid, doi ni el orden de las filas; guardar siempre como .xlsx")))
  openxlsx::setColWidths(wb, "instrucciones", cols = 1:2, widths = c(14, 110))
  openxlsx::saveWorkbook(wb, f)
  invisible(f)
}
walk(1:2, \(a) write_annotator_xlsx(template, a))
# Agreement: irr::kappa2() on direction between the two files (joined by pmid/doi); disagreements resolved by
# consensus -> direction_final.csv (read by 06)

save_versioned(corpus, "corpus_flagged")
