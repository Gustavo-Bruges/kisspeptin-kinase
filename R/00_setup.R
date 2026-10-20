## =============================================================================
## 00_setup.R · Paquetes, rutas, parámetros globales
## Ejecutar al inicio de cada sesión: source("R/00_setup.R")
## =============================================================================

## ---- Paquetes ----
cran_pkgs <- c(
  "dplyr", "tidyr", "stringr", "stringi", "purrr", "readr", "tibble",
  "httr2", "openalexR", "bibliometrix", "igraph", "Matrix", "backbone",
  "babelgene", "rentrez", "irr", "mclust", "ggplot2", "patchwork", "here"
)
bioc_pkgs <- c(
  "GEOquery", "limma", "edgeR", "decoupleR", "OmnipathR", "CARNIVAL",
  "AnnotationDbi", "org.Hs.eg.db"
)

install_if_missing <- function() {
  missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
  if (length(missing_cran)) install.packages(missing_cran)
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
  if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE)
}
# install_if_missing()   # run once; then renv::snapshot() to freeze versions

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(stringi)
  library(purrr); library(readr); library(tibble); library(here)
})

## ---- Rutas ----
paths <- list(
  raw_wos      = here("data", "raw", "wos"),
  raw_scopus   = here("data", "raw", "scopus"),
  raw_oa       = here("data", "raw", "openalex"),
  raw_pubtator = here("data", "raw", "pubtator"),
  raw_geo      = here("data", "raw", "geo"),
  curated      = here("data", "raw", "curated"),
  processed    = here("data", "processed"),
  annotation   = here("data", "annotation"),
  figures      = here("output", "figures"),
  tables       = here("output", "tables"),
  supp         = here("output", "supplementary")
)

## ---- Parámetros preregistrados (no modificar tras el registro en OSF) ----
params <- list(
  seed              = 20261014,
  query_version     = "1.1",
  main_base         = "scopus",              # "wos" or "scopus"; single switch read by 03 and 04 (docs/decisions.md, 2026-10-18)
  years             = c(1996, 2025),         # KiSS-1 description to last complete year
  community_min_n   = 20,                    # clusters smaller than this are labelled "other"
  onc_lift_cut      = 2,                     # oncology vocabulary >= 2x the corpus rate -> "oncology" community
  neuro_lift_floor  = 0.5,                   # neuroendocrine vocabulary below half the corpus rate -> "other"
  kin_comm_min_docs = 3,                     # kinase needs >= 3 papers in labelled communities to get a community (06)
  group_min_shared  = 1,                     # shared authors linking two papers into the same author group (06)
  min_support       = 3,                     # minimum co-occurrence count per gene pair
  fdr_alpha         = 0.05,
  lit_high_quantile = 0.75,                  # lit_score cut for "high literature support"
  lit_min_docs      = 3,                     # and at least 3 corpus papers (lit_score on 1-2 papers is noise)
  data_high_frac    = 0.5,                   # fraction of perturbed series for "high data support"
  stability_frac    = 0.8,                   # subsampling fraction in CARNIVAL robustness
  stability_keep    = 0.7,                   # node retained if present in >= 70% of runs
  n_stability       = 100,
  carnival_timelimit = 600,                    # seconds per CBC run (CARNIVAL default)
  carnival_threads   = 8,                      # CBC threads
  carnival_solver    = "highs",                # "highs" (open source), "gurobi" (academic licence) or "cbc" (optimum only up to depth 4)
  pkn_max_depth      = 5,                      # max path length perturbed node -> TF in the pruned PKN; canonical
                                               # kinases sit 4 steps from KISS1R: 5 is the minimum that lets them reach a TF
                                               # (4 excluded them by construction; 6: no HiGHS optimum in 600 s)
  n_perm            = 10000,                 # transfer test permutations
  periods           = list(p1 = c(1996, 2008), p2 = c(2009, 2016), p3 = c(2017, 2025))
)
set.seed(params$seed)

## ---- Identificación ante APIs ----
# On Windows "~" is Documents (OneDrive here) in Rgui/RStudio but the user folder under Git Bash:
# if the key is missing, load whichever .Renviron exists. Keys are only read, never written or printed.
if (!nzchar(Sys.getenv("OPENALEX_API_KEY"))) {
  renv_files <- unique(c(path.expand("~/.Renviron"),
                         file.path(Sys.getenv("USERPROFILE"), ".Renviron"),
                         file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Documents", ".Renviron")))
  for (f in renv_files[file.exists(renv_files)]) readRenviron(f)
  if (!nzchar(Sys.getenv("OPENALEX_API_KEY"))) warning("OPENALEX_API_KEY not found in: ", paste(renv_files, collapse = ", "), call. = FALSE)
}
options(
  openalexR.mailto = Sys.getenv("OPENALEX_MAILTO"),  # contact email for the OpenAlex polite pool, also in ~/.Renviron
  openalexR.apikey = Sys.getenv("OPENALEX_API_KEY")  # stored in ~/.Renviron, never in the script
)

## ---- Utilidad: guardar con metadatos de reproducibilidad ----
save_versioned <- function(obj, name, ...) {
  saveRDS(
    list(data = obj, saved_at = Sys.time(), query_version = params$query_version,
         r_version = R.version.string, extra = list(...)),
    file.path(paths$processed, paste0(name, ".rds"))
  )
}
load_versioned <- function(name) readRDS(file.path(paths$processed, paste0(name, ".rds")))$data

## ---- Utilidad: llave de registro (DOI; UT si falta el DOI), compartida por 02 y 04 ----
record_key <- function(doi, ut) coalesce(na_if(doi, ""), str_c("ut:", ut))
