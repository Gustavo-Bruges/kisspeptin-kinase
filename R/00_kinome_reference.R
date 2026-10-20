## =============================================================================
## 00_kinome_reference.R · Referencia del kinoma humano (UniProt pkinfam)
## Insumo:  data/raw/curated/pkinfam_<fecha>.txt (UniProt, docs/pkinfam.txt)
## Salida:  data/raw/curated/kinome_reference.csv (columna symbol + familia, accesión, release)
## Fuente elegida: UniProt pkinfam (versionada; incluye quinasas atípicas) · decisions.md 2026-10-18
## =============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })

pk_file <- sort(list.files(paths$curated, "^pkinfam_.*[.]txt$", full.names = TRUE), decreasing = TRUE)[1]
lines   <- readLines(pk_file, warn = FALSE)
release <- str_squish(str_remove(str_subset(lines, "^Release:")[1], "^Release:"))

## ---- Familia vigente para cada línea: el encabezado más reciente subrayado con "====" ----
is_heading <- c(str_detect(lines[-1], "^=+$"), FALSE) & nzchar(str_trim(lines))
family     <- if_else(is_heading, lines, NA_character_)
family     <- tidyr::fill(tibble(family), family)$family

kin <- tibble(line = lines, family = family) |>
  filter(str_detect(line, "_HUMAN +[(]")) |>
  mutate(gene_pkinfam = str_extract(line, "^[^ ]+"),
         uniprot      = str_match(line, "_HUMAN +[(]([A-Z0-9]+)")[, 2],
         group        = str_remove(family, "^Atypical: ") |> str_remove(" (Ser/Thr |Tyr )?(protein )?kinase family$"),
         atypical     = str_detect(family, "^Atypical")) |>
  dplyr::select(gene_pkinfam, uniprot, family, group, atypical) |>   # AnnotationDbi masks select()
  distinct(uniprot, .keep_all = TRUE)

## ---- Símbolo HGNC vigente desde la accesión (pkinfam puede traer nombres antiguos) ----
hgnc <- AnnotationDbi::select(org.Hs.eg.db, keys = kin$uniprot, keytype = "UNIPROT", columns = "SYMBOL") |>
  filter(!is.na(SYMBOL)) |> group_by(UNIPROT) |>
  filter(n_distinct(SYMBOL) == 1) |> ungroup() |> distinct(uniprot = UNIPROT, symbol = SYMBOL)

kinome <- kin |>
  left_join(hgnc, by = "uniprot") |>
  mutate(symbol_source = if_else(is.na(symbol), "pkinfam", "org.Hs.eg.db"),
         symbol        = coalesce(symbol, gene_pkinfam),
         release       = release, source = basename(pk_file)) |>
  distinct(symbol, .keep_all = TRUE) |>
  relocate(symbol)

renamed <- filter(kinome, symbol != gene_pkinfam)
message(sprintf("kinome: %d human kinases (%s); %d renamed to current HGNC symbol; %d without HGNC match",
                nrow(kinome), release, nrow(renamed), sum(kinome$symbol_source == "pkinfam")))
write_csv(kinome, file.path(paths$curated, "kinome_reference.csv"))
