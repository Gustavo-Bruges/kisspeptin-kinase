## =============================================================================
## 06_cooccurrence.R · Coocurrencia de genes (oración = unidad principal)
## Salidas: lit_kin (soporte en literatura por quinasa), redes por estrato,
##          matriz gen × contexto
## =============================================================================
source(here::here("R", "00_setup.R"))
library(Matrix)
library(igraph)

corpus  <- load_versioned("corpus_flagged")
pt_text <- load_versioned("pubtator_text")
genes_h <- load_versioned("genes_human")

kinome     <- read_csv(file.path(paths$curated, "kinome_reference.csv"))$symbol  # fixed reference, version in decisions.md
kiss_axis  <- c("KISS1", "KISS1R")

## ---- A. Segmentación en oraciones y asignación de menciones ----
sentences <- pt_text |>
  mutate(spans = stri_locate_all_boundaries(full, type = "sentence")) |>   # ICU boundaries, 1-based
  select(pmid, spans) |>
  mutate(spans = map(spans, ~ as_tibble(.x) |> mutate(sid = row_number()))) |>
  unnest(spans) |>
  rename(s_start = start, s_end = end)

genes_sent <- genes_h |>
  mutate(pos = start + 1L) |>                                 # PubTator offsets are 0-based
  inner_join(sentences, join_by(pmid, between(pos, s_start, s_end))) |>
  mutate(unit = str_c(pmid, "_", sid)) |>
  distinct(pmid, unit, symbol)

genes_abs <- genes_h |> distinct(pmid, symbol) |> mutate(unit = pmid)   # abstract-level units (sensitivity)

## ---- B. Coocurrencia con NPMI, association strength e hipergeométrica ----
cooc <- function(unit_gene, min_support = params$min_support) {
  N    <- n_distinct(unit_gene$unit)
  freq <- count(unit_gene, symbol, name = "n_g")
  unit_gene |>
    select(unit, symbol) |>
    inner_join(select(unit_gene, unit, symbol), by = "unit", relationship = "many-to-many") |>
    filter(symbol.x < symbol.y) |>
    count(symbol.x, symbol.y, name = "n_xy") |>
    filter(n_xy >= min_support) |>
    left_join(freq, by = c("symbol.x" = "symbol")) |> rename(n_x = n_g) |>
    left_join(freq, by = c("symbol.y" = "symbol")) |> rename(n_y = n_g) |>
    mutate(
      npmi  = log((n_xy / N) / ((n_x / N) * (n_y / N))) / -log(n_xy / N),
      as_ve = n_xy / (n_x * n_y),                            # Van Eck association strength (VOSviewer comparability)
      p     = phyper(n_xy - 1, n_x, N - n_x, n_y, lower.tail = FALSE),
      fdr   = p.adjust(p, method = "BH"),
      N_units = N
    )
}

## ---- C. Backbone bipartito SDSM (grados preservados) ----
to_incidence <- function(unit_gene) {
  ug <- distinct(unit_gene, unit, symbol)
  sparseMatrix(i = as.integer(factor(ug$unit)), j = as.integer(factor(ug$symbol)), x = 1,
               dimnames = list(levels(factor(ug$unit)), levels(factor(ug$symbol))))
}
# backbone >= 3.0: sdsm() is deprecated and lost its `class` argument; backbone_from_projection(model = "sdsm")
# takes a sparse genes x units matrix (it projects onto ROWS) and returns the gene x gene backbone adjacency.
# Genes seen in fewer than min_support units cannot reach min_support co-occurrences and are dropped first.
backbone_edges <- function(unit_gene) {
  keep <- unit_gene |> distinct(unit, symbol) |> count(symbol) |> filter(n >= params$min_support) |> pull(symbol)
  B  <- Matrix::t(to_incidence(filter(unit_gene, symbol %in% keep)))
  bb <- backbone::backbone_from_projection(B, alpha = params$fdr_alpha, mtc = "fdr", model = "sdsm")
  e  <- summary(as(as(bb, "generalMatrix"), "TsparseMatrix"))
  tibble(symbol.x = rownames(B)[e$i], symbol.y = rownames(B)[e$j]) |> filter(symbol.x < symbol.y)
}

## ---- D. Red principal (oración, todo el corpus) ----
pairs_sent <- cooc(genes_sent)
bb_sent    <- backbone_edges(genes_sent)

edges_main <- pairs_sent |>
  semi_join(bind_rows(bb_sent, rename(bb_sent, symbol.x = symbol.y, symbol.y = symbol.x)),
            by = c("symbol.x", "symbol.y")) |>
  filter(fdr < params$fdr_alpha, npmi > 0,
         symbol.x %in% kinome | symbol.y %in% kinome)        # edges touching at least one kinase
write_csv(edges_main, file.path(paths$tables, "edges_sentence_main.csv"))

## Sensibilidad: unidad = resumen; y sin revisiones (enumeran genes sin relacionarlos)
pairs_abs       <- cooc(genes_abs)
pairs_noreview  <- cooc(genes_sent |> filter(pmid %in% corpus$pmid[!corpus$is_review]))
write_csv(pairs_abs,      file.path(paths$supp, "S4_pairs_abstract_level.csv"))
write_csv(pairs_noreview, file.path(paths$supp, "S4_pairs_no_reviews.csv"))

## ---- E. Soporte en literatura por quinasa (lit_kin) ----
# Background = NCBI gene2pubmed (decisions.md 2026-10-18). It is gene-normalized like PubTator (ERK2 and MAPK1
# share one GeneID), so aliases no longer inflate lit_score. Counts pool the human gene and its NCBI orthologs,
# because 05 maps mouse/rat/sheep mentions to the human symbol. Universe on each side: PMIDs with any gene
# (gene2pubmed for PubMed, PubTator for the corpus). gene2pubmed mixes curated and text-mined links, so it is
# not PubTator itself: a residual difference in sensitivity remains and is stated as a limitation.
lit_background_method <- "gene2pubmed_with_orthologs"
kin_ref    <- read_csv(file.path(paths$curated, "kinome_reference.csv"), show_col_types = FALSE)
kin_entrez <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = kin_ref$symbol, keytype = "SYMBOL",
                                    columns = "ENTREZID") |>
  filter(!is.na(ENTREZID)) |> transmute(symbol = SYMBOL, human_entrez = ENTREZID)
ncbi_pairs <- readRDS(file.path(paths$curated, "gene_orthologs_human_pairs.rds"))$pairs
kin_genes  <- bind_rows(transmute(kin_entrez, symbol, gene = human_entrez),
                        inner_join(kin_entrez, ncbi_pairs, by = "human_entrez") |> transmute(symbol, gene = entrez)) |>
  mutate(gene = as.integer(gene)) |> distinct()

g2p_raw   <- sort(list.files(paths$curated, "^gene2pubmed_.*[.]gz$", full.names = TRUE), decreasing = TRUE)[1]
g2p_cache <- file.path(paths$curated, "gene2pubmed_kinome_counts.rds")
if (!file.exists(g2p_cache) || !identical(readRDS(g2p_cache)$key, c(basename(g2p_raw), unique(kin_ref$release)))) {
  g2p <- data.table::fread(g2p_raw, select = c(2L, 3L), col.names = c("gene", "pmid"))   # 84M rows, ~8 s
  N_bg <- data.table::uniqueN(g2p$pmid)                                                   # PubMed papers with any gene
  bg   <- as_tibble(g2p[g2p$gene %in% kin_genes$gene, ]) |> inner_join(kin_genes, by = "gene",
                                                                         relationship = "many-to-many") |>
    group_by(symbol) |> summarise(pubmed_n = n_distinct(pmid), .groups = "drop")
  saveRDS(list(bg = bg, N_bg = N_bg, key = c(basename(g2p_raw), unique(kin_ref$release))), g2p_cache)
  rm(g2p)
}
g2p_obj  <- readRDS(g2p_cache)
N_pubmed <- g2p_obj$N_bg
bg       <- g2p_obj$bg

N_corpus <- n_distinct(genes_abs$pmid)                          # corpus papers with any PubTator gene
kiss_sent_units <- genes_sent |> filter(symbol %in% kiss_axis) |> distinct(unit)

## Community of each kinase = the one where it is most over-represented (lift), not the most frequent:
## 85 % of the corpus is neuroendocrine, so "most frequent" labels almost every kinase neuroendocrine
labelled   <- c("oncology", "neuroendocrine", "molecular")
comm_share <- corpus |> filter(pmid %in% genes_abs$pmid, community %in% labelled) |>
  count(community, name = "n_comm") |> mutate(p_comm = n_comm / sum(n_comm))
kin_comm <- genes_abs |> filter(symbol %in% kinome) |> distinct(symbol, pmid) |>
  inner_join(select(corpus, pmid, community), by = "pmid") |> filter(community %in% labelled) |>
  count(symbol, community) |> group_by(symbol) |> mutate(n_labelled = sum(n), frac = n / n_labelled) |> ungroup() |>
  inner_join(comm_share, by = "community") |> mutate(lift = frac / p_comm)
write_csv(kin_comm, file.path(paths$supp, "S4_kinase_community_lift.csv"))
kin_comm_best <- kin_comm |> group_by(symbol) |> slice_max(lift, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(symbol, community = if_else(n_labelled >= params$kin_comm_min_docs, community, "insufficient"),
            community_lift = lift)

lit_kin <- genes_abs |>
  filter(symbol %in% kinome) |>
  left_join(select(corpus, pmid, year), by = "pmid") |>
  group_by(symbol) |>
  summarise(n_docs = n_distinct(pmid), first_year = min(year, na.rm = TRUE), .groups = "drop") |>
  left_join(kin_comm_best, by = "symbol") |>
  mutate(community = coalesce(community, "insufficient")) |>
  left_join(genes_sent |> filter(symbol %in% kinome) |> semi_join(kiss_sent_units, by = "unit") |>
              count(symbol, name = "n_sent_kiss"), by = "symbol") |>
  left_join(bg, by = "symbol") |>
  mutate(
    n_sent_kiss = coalesce(n_sent_kiss, 0L),
    pubmed_n    = coalesce(pubmed_n, 0L),
    lit_score   = log2(((n_docs + 1) / N_corpus) / ((pubmed_n + 1) / N_pubmed))  # enrichment vs gene2pubmed background
  )
save_versioned(lit_kin, "lit_kin", N_pubmed = N_pubmed, N_corpus = N_corpus, background = lit_background_method,
               gene2pubmed = basename(g2p_raw))
write_csv(lit_kin, file.path(paths$tables, "lit_kin.csv"))

## ---- F. Anotación con recursos curados (capa, no filtro) ----
# Optional until downloaded (CGC needs a COSMIC licence): missing layers come out as NA, never as FALSE
read_optional <- function(f, reader) if (file.exists(f)) reader(f, show_col_types = FALSE) else {
  warning("annotation layer missing, columns left NA: ", basename(f), call. = FALSE); NULL }
cgc        <- read_optional(file.path(paths$curated, "cancer_gene_census.csv"), read_csv)   # COSMIC release recorded
cancermine <- read_optional(file.path(paths$curated, "cancermine_collated.tsv"), read_tsv)

gene_annot <- tibble(symbol = unique(genes_sent$symbol)) |>
  mutate(
    is_kinase     = symbol %in% kinome,
    in_cgc        = if (is.null(cgc)) NA else symbol %in% cgc$`Gene Symbol`,
    cgc_role      = if (is.null(cgc)) NA_character_ else cgc$`Role in Cancer`[match(symbol, cgc$`Gene Symbol`)],  # multi-role kept as listed
    in_cancermine = if (is.null(cancermine)) NA else symbol %in% cancermine$gene_normalized
  )
write_csv(gene_annot, file.path(paths$tables, "gene_annotation.csv"))

## ---- G. Estratos y matriz gen × contexto ----
dir_file <- file.path(paths$annotation, "direction_final.csv")   # produced after double annotation
dir_final <- if (file.exists(dir_file)) read_csv(dir_file, col_types = cols(pmid = "c")) else {
  warning("direction_final.csv not available yet: onc_suppr stratum and breast_promo table are empty", call. = FALSE)
  tibble(pmid = character(), doi = character(), year = integer(), title = character(), direction = character(),
         model = character(), tnbc = logical(), breast = logical(), er_neg = logical())
}

strata <- list(
  all          = corpus$pmid,
  neuroendo    = corpus$pmid[corpus$community == "neuroendocrine"],
  oncology     = corpus$pmid[corpus$community == "oncology"],
  molecular    = corpus$pmid[corpus$community == "molecular"],
  onc_suppr    = dir_final$pmid[dir_final$direction == "suppressive"],
  p1           = corpus$pmid[between(corpus$year, params$periods$p1[1], params$periods$p1[2])],
  p2           = corpus$pmid[between(corpus$year, params$periods$p2[1], params$periods$p2[2])],
  p3           = corpus$pmid[between(corpus$year, params$periods$p3[1], params$periods$p3[2])]
)
# breast_promo (promotive & (TNBC | breast ER-negative)) is NOT a network stratum: ~4-5 papers, mostly one
# laboratory (decisions.md 2026-10-18). It is reported descriptively in section I.
breast_promo <- dir_final$pmid[dir_final$direction == "promotive" &
                                 (dir_final$tnbc | (dir_final$breast & dir_final$er_neg))]   # ER-negative endometrial etc. excluded

# Gene x context matrix = prevalence of each kinase per stratum and its lift against the whole corpus.
# NPMI with KISS1/KISS1R was dropped (decisions.md 2026-10-18): the axis is in 54 % of gene sentences and 71 %
# of papers of a corpus selected for it, so co-occurring "with KISS1" carries almost no information.
gene_docs <- genes_abs |> distinct(pmid, symbol)
prev_all  <- gene_docs |> filter(symbol %in% kinome) |> count(symbol, name = "n_all") |>
  mutate(prev_all = n_all / n_distinct(gene_docs$pmid))
gene_ctx_long <- imap_dfr(strata, \(ids, s) {
  docs <- intersect(unique(na.omit(ids)), unique(gene_docs$pmid))        # stratum papers with any gene
  gene_docs |> filter(pmid %in% docs, symbol %in% kinome) |> count(symbol, name = "n_docs") |>
    mutate(stratum = s, N_stratum = length(docs), prevalence = n_docs / N_stratum)
}) |>
  left_join(prev_all, by = "symbol") |>
  mutate(lift = prevalence / prev_all) |>
  select(symbol, stratum, n_docs, N_stratum, prevalence, lift)
write_csv(gene_ctx_long, file.path(paths$supp, "S4_kinase_prevalence_by_stratum.csv"))
gene_ctx <- gene_ctx_long |>
  filter(n_docs >= params$min_support) |>                     # lifts on 1-2 papers are not reported
  select(symbol, stratum, lift) |>
  pivot_wider(names_from = stratum, values_from = lift)
write_csv(gene_ctx, file.path(paths$tables, "gene_by_context_lift.csv"))

## ---- H. Redes por estrato para figuras ----
nets <- imap(strata[c("all", "neuroendo", "oncology", "onc_suppr")], \(ids, s) {
  e <- cooc(filter(genes_sent, pmid %in% na.omit(ids))) |>
    filter(fdr < params$fdr_alpha, npmi > 0, symbol.x %in% kinome | symbol.y %in% kinome)
  if (nrow(e) == 0) { warning("stratum ", s, ": no significant kinase edges, network skipped", call. = FALSE); return(NULL) }
  g <- graph_from_data_frame(select(e, symbol.x, symbol.y, weight = npmi), directed = FALSE)
  V(g)$is_kinase <- V(g)$name %in% kinome
  set.seed(params$seed)
  V(g)$cluster <- membership(cluster_louvain(g, weights = E(g)$weight))
  g
})
saveRDS(nets, file.path(paths$processed, "networks_by_stratum.rds"))

## ---- I. Grupos de autores por estrato y vista descriptiva breast_promo ----
# A small stratum can be a single laboratory. Author groups = connected components of the paper graph where
# two papers are linked if they share >= params$group_min_shared authors (bibliometrix AU: "SURNAME I").
# Homonyms ("WANG Y") can merge distinct groups, so n_groups is a conservative (lower) count.
M_main  <- load_versioned(paste0("M_", params$main_base))
authors <- M_main |> transmute(key = record_key(DI, UT), AU) |>
  left_join(select(corpus, key, pmid), by = "key") |>        # key (DOI, else UT): records without DOI carry "" and would cross-match
  filter(!is.na(pmid)) |> separate_rows(AU, sep = ";") |> mutate(AU = str_squish(AU)) |>
  filter(AU != "") |> distinct(pmid, AU)

author_groups <- function(ids) {
  a <- filter(authors, pmid %in% ids)
  if (n_distinct(a$pmid) == 0) return(tibble(n_papers = 0L, n_groups = 0L, largest_group_share = NA_real_))
  links <- a |> inner_join(a, by = "AU", relationship = "many-to-many") |> filter(pmid.x < pmid.y) |>
    count(pmid.x, pmid.y) |> filter(n >= params$group_min_shared)
  g  <- graph_from_data_frame(select(links, pmid.x, pmid.y), directed = FALSE,
                              vertices = tibble(name = unique(a$pmid)))
  cc <- components(g)
  tibble(n_papers = vcount(g), n_groups = cc$no, largest_group_share = max(cc$csize) / vcount(g))
}
groups_by_stratum <- imap_dfr(c(strata, list(breast_promo = breast_promo)),
                              \(ids, s) author_groups(na.omit(ids)) |> mutate(stratum = s, .before = 1))
write_csv(groups_by_stratum, file.path(paths$tables, "author_groups_by_stratum.csv"))

kin_by_paper <- genes_abs |> filter(symbol %in% kinome) |> group_by(pmid) |>
  summarise(kinases = paste(sort(unique(symbol)), collapse = "; "), .groups = "drop")
breast_promo_tbl <- dir_final |> filter(pmid %in% breast_promo) |>
  select(pmid, doi, year, title, direction, model, any_of("notes")) |>
  left_join(kin_by_paper, by = "pmid") |>
  left_join(authors |> group_by(pmid) |> summarise(authors = paste(AU, collapse = "; "), .groups = "drop"), by = "pmid") |>
  arrange(year)
write_csv(breast_promo_tbl, file.path(paths$tables, "breast_promo_descriptive.csv"))
