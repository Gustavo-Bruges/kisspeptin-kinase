## =============================================================================
## 01_search_openalex.R · Corpus OpenAlex (análisis de sensibilidad entre bases)
## Ecuación v1.1 · docs/search_equations_v1.1.txt
## =============================================================================
source(here::here("R", "00_setup.R"))
library(openalexR)

stopifnot(nzchar(Sys.getenv("OPENALEX_API_KEY")))  # cursor paging requires a free API key

## ---- Ecuación C (sin comodines: OpenAlex no los expande) ----
kp_query <- paste(
  '(kisspeptin OR kisspeptins OR "kisspeptin-10" OR "kisspeptin-54" OR metastin',
  'OR KISS1 OR KISS1R OR "KiSS-1" OR "KiSS 1 receptor" OR kiss2 OR kissr2 OR kissr3',
  'OR GPR54 OR "GPR 54" OR AXOR12 OR hOT7T175)'
)

oa_args <- list(
  entity = "works",
  title_and_abstract.search = kp_query,
  from_publication_date = sprintf("%d-01-01", params$years[1]),
  to_publication_date   = sprintf("%d-12-31", params$years[2]),
  type = c("article", "review")
)

## ---- Descarga con caché (re-ejecutar no vuelve a consultar la API) ----
# Raw download kept in data/raw/openalex; delete the file to force a fresh extraction (new extraction date)
raw_file <- file.path(paths$raw_oa, sprintf("openalex_raw_v%s.rds", params$query_version))

if (!file.exists(raw_file)) {
  print(do.call(oa_query, oa_args))                          # paste in browser to inspect parsing
  n_expected <- do.call(oa_fetch, c(oa_args, count_only = TRUE))$count
  message("OpenAlex expected count: ", n_expected)
  raw <- do.call(oa_fetch, c(oa_args, abstract = TRUE, verbose = TRUE))
  saveRDS(list(data = raw, fetched_at = Sys.time(), n_expected = n_expected,
               openalexR = as.character(packageVersion("openalexR"))), raw_file)
}
raw_obj    <- readRDS(raw_file)
corpus_oa  <- raw_obj$data
n_expected <- raw_obj$n_expected

## ---- Armonización de columnas (cambian entre versiones de openalexR) ----
# openalexR >= 2.0 already returns title and abstract (plus display_name); renaming blindly duplicates `title`
if (!"title" %in% names(corpus_oa))    corpus_oa <- rename(corpus_oa, title = display_name)
if (!"abstract" %in% names(corpus_oa)) corpus_oa <- rename(corpus_oa, abstract = ab)

n_fetched <- nrow(corpus_oa)
corpus_oa <- corpus_oa |>
  distinct(id, .keep_all = TRUE) |>                          # cursor paging can return the same work twice
  mutate(
    doi    = str_to_lower(str_remove(doi, "^https?://doi.org/")),  # bare lowercase DOI as cross-database key
    source = "openalex"
  )

## ---- Control de calidad ----
n_refs <- map_int(corpus_oa$referenced_works, \(x) sum(!is.na(x)))   # empty list entries come back as NA, not length 0
qc_oa <- corpus_oa |>
  summarise(
    n_expected  = n_expected,
    n_fetched   = n_fetched,
    n           = n(),                                       # after removing duplicated work ids
    pct_no_abs  = mean(is.na(abstract) | abstract == "") * 100,
    pct_no_doi  = mean(is.na(doi)) * 100,
    pct_no_refs = mean(n_refs == 0) * 100
  )
print(qc_oa)
write_csv(qc_oa, file.path(paths$supp, "S1_qc_openalex.csv"))

save_versioned(corpus_oa, "corpus_openalex",
               query = kp_query, n_expected = n_expected, fetched_at = raw_obj$fetched_at,
               openalexR = raw_obj$openalexR)
