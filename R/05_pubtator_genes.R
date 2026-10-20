## =============================================================================
## 05_pubtator_genes.R · Anotación de genes (PubTator3) y símbolo humano oficial
## Unidad de anotación: título + resumen
## =============================================================================
source(here::here("R", "00_setup.R"))
library(httr2)

corpus <- load_versioned("corpus_flagged")
pmids  <- unique(na.omit(corpus$pmid))   # records without PMID cannot be annotated: loss reported in 02

## ---- Descarga con caché (re-ejecutar no vuelve a consultar la API) ----
cache_file <- file.path(paths$raw_pubtator, "pubtator_raw_v1.1.txt")

fetch_pubtator <- function(ids) {
  # endpoint and batch limit per PubTator3 docs at time of writing; verify before batch use
  request("https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/pubtator") |>
    req_url_query(pmids = paste(ids, collapse = ",")) |>
    req_retry(max_tries = 3) |>
    req_throttle(rate = 3 / 1) |>                      # NCBI rate limit (~3 requests/s)
    req_perform() |>
    resp_body_string()
}

if (!file.exists(cache_file)) {
  batches <- split(pmids, ceiling(seq_along(pmids) / 100))
  raw_txt <- map_chr(batches, fetch_pubtator, .progress = TRUE)
  writeLines(raw_txt, cache_file)
}
raw_lines <- readLines(cache_file)

## ---- Texto (título y resumen) para segmentar oraciones en 06 ----
pt_text <- tibble(line = raw_lines) |>
  filter(str_detect(line, "^\\d+\\|[ta]\\|")) |>
  extract(line, c("pmid", "part", "text"), "^(\\d+)\\|([ta])\\|(.*)$") |>
  pivot_wider(names_from = part, values_from = text) |>
  mutate(full = str_c(t, " ", coalesce(a, "")))     # PubTator offsets index title + single space + abstract

## ---- Anotaciones ----
ann <- tibble(line = raw_lines) |>
  filter(str_count(line, "\t") >= 5) |>
  separate(line, into = c("pmid", "start", "end", "mention", "type", "id"),
           sep = "\t", extra = "drop", convert = TRUE) |>
  mutate(pmid = as.character(pmid))                  # character key for joins with corpus and sentences

genes <- ann |>
  filter(type == "Gene", !is.na(id), id != "") |>
  separate_rows(id, sep = ";") |>                    # one mention may carry several Gene IDs
  mutate(entrez = as.character(id))

## ---- Símbolo humano oficial: humanos directos + ortólogos ----
library(AnnotationDbi)
library(org.Hs.eg.db)
library(babelgene)                                   # orthologs() accepts Entrez IDs (verified with babelgene 22.9)

human_map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(genes$entrez),
                                   keytype = "ENTREZID", columns = "SYMBOL") |>
  filter(!is.na(SYMBOL)) |>                          # non-human Entrez IDs return NA here
  transmute(entrez = ENTREZID, human_symbol = SYMBOL)

nonhuman_ids <- setdiff(unique(genes$entrez), human_map$entrez)

# Every species bundled with babelgene (19 species). Sheep, goat, hamster, medaka, tilapia... are not among
# them, and babelgene lacks KISS1 even for some it covers (pig, dog, cat).
orth_species <- setdiff(babelgene::species()$scientific_name, "Homo sapiens")
orth_map <- map_dfr(orth_species, \(sp) {
  tryCatch(orthologs(genes = nonhuman_ids, species = sp, human = FALSE),
           error = function(e) NULL)                 # species without matches returns nothing
}) |>
  transmute(entrez = as.character(entrez), human_symbol) |>   # babelgene columns: symbol, entrez, human_symbol, human_entrez
  distinct()

## ---- Tercera fuente: ortología de NCBI Gene (gene_orthologs, todas las especies) ----
# First pass left ~10 % of mentions unmapped, half of them sheep/goat KISS1, GNRH1, TAC3, KISS1R: the ruminant
# community would lose its own axis genes. NCBI pairs every species with a human gene where one exists.
# The raw file (data/raw/curated/gene_orthologs_<date>.gz) is reduced once to human pairs and cached.
ncbi_raw   <- sort(list.files(paths$curated, "^gene_orthologs_.*[.]gz$", full.names = TRUE), decreasing = TRUE)[1]
ncbi_cache <- file.path(paths$curated, "gene_orthologs_human_pairs.rds")
if (!file.exists(ncbi_cache)) {
  go <- read_tsv(ncbi_raw, col_names = c("tax_id", "gene", "relationship", "other_tax_id", "other_gene"),
                 skip = 1, col_types = cols_only(tax_id = "i", gene = "c", other_tax_id = "i", other_gene = "c"))
  pairs <- bind_rows(
    go |> filter(other_tax_id == 9606) |> transmute(entrez = gene, tax_id, human_entrez = other_gene),
    go |> filter(tax_id == 9606)       |> transmute(entrez = other_gene, tax_id = other_tax_id, human_entrez = gene)
  ) |> distinct()
  saveRDS(list(pairs = pairs, genes_with_orthologs = unique(c(go$gene, go$other_gene)),
               source = basename(ncbi_raw)), ncbi_cache)
  rm(go)
}
ncbi <- readRDS(ncbi_cache)

still_unmapped <- setdiff(nonhuman_ids, orth_map$entrez)
ncbi_hits <- ncbi$pairs |> filter(entrez %in% still_unmapped)
ncbi_sym  <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ncbi_hits$human_entrez),
                                   keytype = "ENTREZID", columns = "SYMBOL") |>
  transmute(human_entrez = ENTREZID, human_symbol = SYMBOL)
ncbi_map <- ncbi_hits |> inner_join(ncbi_sym, by = "human_entrez") |> filter(!is.na(human_symbol)) |>
  group_by(entrez) |> filter(n_distinct(human_symbol) == 1) |> ungroup() |>   # one-to-many orthology left unmapped
  distinct(entrez, human_symbol)

map_all <- bind_rows(mutate(human_map, map_source = "human"),
                     mutate(orth_map,  map_source = "babelgene"),
                     mutate(ncbi_map,  map_source = "ncbi_gene_orthologs")) |>
  distinct(entrez, .keep_all = TRUE)

genes_h <- genes |>
  left_join(map_all, by = "entrez") |>
  mutate(symbol = human_symbol)                      # official HGNC symbol only; raw mention never becomes a node

## ---- Auditoría S3: especie de cada gen sin mapear y pérdida biológica vs hueco del recurso ----
# no_human_ortholog = TRUE: NCBI lists orthologs for the gene but none in human (teleost kiss2, gnrh3: genes
# lost in the mammalian lineage). FALSE: absent from NCBI orthology altogether (resource gap or obsolete ID).
unmapped_ids <- unique(genes_h$entrez[is.na(genes_h$symbol)])
taxa_cache   <- file.path(paths$raw_pubtator, "unmapped_gene_taxa_v1.1.csv")
if (!file.exists(taxa_cache)) {
  taxa <- map_dfr(split(unmapped_ids, ceiling(seq_along(unmapped_ids) / 200)), \(ids) {
    Sys.sleep(if (nzchar(Sys.getenv("ENTREZ_KEY"))) 0.11 else 0.35)
    r <- rentrez::entrez_summary(db = "gene", id = ids, always_return_list = TRUE)
    tibble(entrez = names(r), organism = map_chr(r, \(x) x$organism$scientificname %||% NA_character_),
           ncbi_symbol = map_chr(r, \(x) x$name %||% NA_character_))
  })
  write_csv(taxa, taxa_cache)
}
taxa <- read_csv(taxa_cache, col_types = cols(.default = "c"))

unmapped <- genes_h |> filter(is.na(symbol)) |> count(entrez, mention, sort = TRUE) |>
  left_join(taxa, by = "entrez") |>
  mutate(no_human_ortholog = entrez %in% ncbi$genes_with_orthologs)
write_csv(unmapped, file.path(paths$supp, "S3_genes_unmapped_audit.csv"))
write_csv(unmapped |> group_by(organism) |>
            summarise(mentions = sum(n), no_human_ortholog = sum(n[no_human_ortholog]), .groups = "drop") |>
            arrange(desc(mentions)),
          file.path(paths$supp, "S3_genes_unmapped_by_species.csv"))
write_csv(count(genes_h, map_source, name = "mentions"), file.path(paths$supp, "S3_gene_mapping_sources.csv"))
message(sprintf("Unmapped gene mentions: %.1f%% (%.1f%% without any human ortholog)",
                mean(is.na(genes_h$symbol)) * 100,
                sum(unmapped$n[unmapped$no_human_ortholog]) / nrow(genes_h) * 100))

genes_h <- filter(genes_h, !is.na(symbol))

save_versioned(pt_text, "pubtator_text")
save_versioned(genes_h, "genes_human")
