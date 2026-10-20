## =============================================================================
## 02_import_wos_scopus.R · Importación, depuración y llave DOI→PMID
## Insumos: data/raw/wos/*.txt (plaintext, full record + cited refs)
##          data/raw/scopus/*.csv (o *.bib, con referencias)
## Bases analizadas POR SEPARADO; nunca se fusionan para co-citación
## =============================================================================
source(here::here("R", "00_setup.R"))
library(bibliometrix)
library(rentrez)

## ---- Importación por lotes de exporte ----
read_base <- function(dir, dbsource, format) {
  files <- list.files(dir, pattern = if (format == "csv") "\\.csv$" else "\\.(txt|bib)$",
                      full.names = TRUE)
  stopifnot(length(files) > 0)
  convert2df(files, dbsource = dbsource, format = format)  # bibliometrix field tags: DI, PY, DT, CR, PM, UT
}

M_wos    <- read_base(paths$raw_wos,    "wos",    "plaintext")
M_scopus <- read_base(paths$raw_scopus, "scopus", "csv")       # switch to "bibtex" if exported as .bib
# bibliometrix 5.5 keeps Scopus "PubMed ID" as PubMed.ID instead of mapping it to PM
if (!"PM" %in% names(M_scopus) && "PubMed.ID" %in% names(M_scopus)) M_scopus <- rename(M_scopus, PM = PubMed.ID)

## ---- Control: exportes con referencias citadas (convert2df no avisa si faltan; 03 fallaría sin explicación) ----
check_refs <- function(M, base, min_frac = 0.9) {
  frac <- mean(!is.na(M$CR) & M$CR != "")
  if (frac < min_frac)
    stop(sprintf("%s: only %.1f%% of records have cited references (CR). Re-export including References.",
                 base, frac * 100), call. = FALSE)
  message(sprintf("%s: %.1f%% of records with cited references", base, frac * 100))
}
check_refs(M_wos, "wos")
check_refs(M_scopus, "scopus")

## ---- Depuración (misma regla en ambas bases) ----
clean_base <- function(M, base) {
  n0 <- nrow(M)
  M <- M |>
    mutate(DI = str_to_lower(str_trim(DI)), base = base) |>
    filter(between(PY, params$years[1], params$years[2])) |>
    filter(str_detect(str_to_lower(DT), "article|review")) |>     # erratum, editorial, meeting abstract excluded
    filter(!str_detect(str_to_lower(DT), "retract|correction|erratum"))
  M_doi   <- M |> filter(!is.na(DI) & DI != "") |> distinct(DI, .keep_all = TRUE)
  M_nodoi <- M |> filter(is.na(DI) | DI == "") |> distinct(UT, .keep_all = TRUE)
  M <- bind_rows(M_doi, M_nodoi)
  attr(M, "flow") <- tibble(base = base, n_raw = n0, n_clean = nrow(M),
                            n_no_doi = nrow(M_nodoi))           # counts for PRISMA-style flow figure
  M
}

M_wos    <- clean_base(M_wos, "wos")
M_scopus <- clean_base(M_scopus, "scopus")

## ---- Exclusión de registros ajenos al tema (lista auditable, confirmada a mano) ----
# Homonyms captured by the truncations (KISSR galaxy catalog, kissra bread, KISS2 FSM format...).
# Only rows with decision == "exclude" AND a reviewer are applied; unconfirmed proposals are reported, not applied.
excl_file <- file.path(paths$annotation, "offtopic_exclusions.csv")
excl <- if (file.exists(excl_file)) read_csv(excl_file, col_types = cols(.default = "c")) else
  tibble(base = character(), id = character(), decision = character(), reviewer = character())
n_pending <- sum(excl$decision == "exclude" & is.na(excl$reviewer), na.rm = TRUE)
if (n_pending > 0) warning(n_pending, " proposed exclusions without reviewer: not applied", call. = FALSE)
excl_ok <- filter(excl, decision == "exclude", !is.na(reviewer))

drop_offtopic <- function(M, base) {
  flow <- attr(M, "flow")
  ids  <- excl_ok$id[excl_ok$base == base]      # computed outside filter(): there `base` would mask to M's column
  M <- filter(M, !UT %in% ids)
  attr(M, "flow") <- mutate(flow, n_offtopic = n_clean - nrow(M), n_clean = nrow(M),
                            n_no_doi = sum(is.na(M$DI) | M$DI == ""))   # flow: raw -> clean -> off-topic removed
  M
}
M_wos    <- drop_offtopic(M_wos, "wos")
M_scopus <- drop_offtopic(M_scopus, "scopus")

flow <- bind_rows(attr(M_wos, "flow"), attr(M_scopus, "flow"))
print(flow)
write_csv(flow, file.path(paths$supp, "S1_flow_counts.csv"))

## ---- Solapamiento entre bases (descriptivo, no fusión) ----
overlap <- tibble(
  only_wos    = length(setdiff(na.omit(M_wos$DI), na.omit(M_scopus$DI))),
  only_scopus = length(setdiff(na.omit(M_scopus$DI), na.omit(M_wos$DI))),
  both        = length(intersect(na.omit(M_wos$DI), na.omit(M_scopus$DI)))
)
write_csv(overlap, file.path(paths$supp, "S1_overlap_doi.csv"))

## ---- Llave DOI → PMID para PubTator3 ----
# key = DOI when present, otherwise "ut:<UT>"; records without DOI keep their native PMID (frequent before 2005)
record_keys <- function(M) {
  pm <- if ("PM" %in% names(M)) as.character(M$PM) else NA_character_   # PM may be absent in Scopus import: filled below
  M |> transmute(key = record_key(DI, UT), doi = na_if(DI, ""),
                 pmid = na_if(str_trim(pm), ""), year = PY, base)
}

records <- bind_rows(record_keys(M_wos), record_keys(M_scopus)) |>
  group_by(key) |>
  summarise(doi = first(doi), pmid = first(na.omit(pmid)), year = first(year),
            bases = paste(sort(unique(base)), collapse = "+"), .groups = "drop")

doi_to_pmid <- function(doi) {
  Sys.sleep(if (nzchar(Sys.getenv("ENTREZ_KEY"))) 0.11 else 0.35)  # NCBI rate limit with/without key
  res <- tryCatch(entrez_search(db = "pubmed", term = paste0(doi, "[doi]"), retmax = 1),
                  error = function(e) NULL)
  if (is.null(res) || length(res$ids) == 0) NA_character_ else res$ids[1]
}

to_resolve <- which(is.na(records$pmid) & !is.na(records$doi))   # no DOI and no PM: cannot be resolved
records$pmid[to_resolve] <- map_chr(records$doi[to_resolve], doi_to_pmid, .progress = TRUE)

pmid_loss <- records |> summarise(n = n(), with_pmid = sum(!is.na(pmid)),
                                  pct_lost = mean(is.na(pmid)) * 100)
print(pmid_loss)                                                  # records without PMID cannot be annotated
write_csv(pmid_loss, file.path(paths$supp, "S1_pmid_coverage.csv"))

writeLines(unique(na.omit(records$pmid)),
           file.path(paths$processed, "corpus_pmids.txt"))     # plain list used by 07 for the OmniPath circularity control

save_versioned(M_wos,    "M_wos")
save_versioned(M_scopus, "M_scopus")
save_versioned(records,  "records_doi_pmid")
