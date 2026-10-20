## =============================================================================
## 03_bibliometrics.R · Comunidades, conectividad, temas, sensibilidad entre bases
## H1: comunidades oncológica y neuroendocrina con baja conectividad
## =============================================================================
source(here::here("R", "00_setup.R"))
library(bibliometrix)
library(igraph)
library(Matrix)

corpus_oa <- load_versioned("corpus_openalex")

## ---- Base principal (params$main_base, fijada en docs/decisions.md antes de correr) ----
M_main <- load_versioned(paste0("M_", params$main_base))

## ---- Comunidades de artículos del corpus por acoplamiento bibliográfico ----
# Bibliographic coupling clusters the corpus papers themselves (each paper gets a label);
# co-citation clusters the cited references (intellectual base) and is reported separately.
louvain_communities <- function(adj, ids) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_vertices(g, degree(g) == 0)                   # isolates cannot be assigned
  set.seed(params$seed)
  cl <- cluster_louvain(g, weights = E(g)$weight)
  tibble(id = V(g)$name, cluster = as.integer(membership(cl)))
}

# Association strength for coupling: S_ij = C_ij / (k_i * k_j), k = number of references of each paper.
# Sparse diagonal scaling keeps the 4000 x 4000 matrix sparse (no dense outer product)
association_strength <- function(B) {
  k <- Matrix::rowSums(B)
  D <- Diagonal(x = 1 / k)
  S <- D %*% tcrossprod(B) %*% D                            # shared references between papers, normalized
  dimnames(S) <- list(rownames(B), rownames(B))
  as(forceSymmetric(S, uplo = "U"), "generalMatrix")        # D C D is symmetric up to rounding; igraph requires exact symmetry
}

# WoS / Scopus: paper x cited-reference incidence from bibliometrix, coupling computed here.
# biblioNetwork(analysis = "coupling") is not used: for Scopus it keeps only rows whose name starts with a
# letter (meant for cited-reference labels), which silently drops every paper named by its DOI ("10.")
coupling_bib <- function(M) {
  M <- M |> filter(!is.na(DI) & DI != "")                  # DOI is the cross-database key; records without DOI excluded from ARI
  M <- M |> filter(!is.na(CR) & CR != "")                  # no references: zero self-coupling gives NaN under association strength
  M$SR <- M$DI                                              # cocMatrix names rows by SR; DOI used instead
  rownames(M) <- M$SR
  B <- cocMatrix(M, Field = "CR", type = "sparse", sep = ";")
  B@x[] <- 1                                                # presence, not count, of each cited reference
  B <- B[, nchar(colnames(B)) > 0, drop = FALSE]
  association_strength(B[Matrix::rowSums(B) > 0, , drop = FALSE])
}
adj_main  <- coupling_bib(M_main)                          # computed once; reused for the connectivity index
comm_main <- louvain_communities(adj_main, M_main$DI) |> rename(doi = id)

# OpenAlex: coupling from referenced_works identifiers (no string disambiguation needed)
coupling_oa <- function(df) {
  long <- df |> select(doi, referenced_works) |> filter(!is.na(doi)) |>
    unnest_longer(referenced_works) |> filter(!is.na(referenced_works)) |> distinct()
  B <- sparseMatrix(i = as.integer(factor(long$doi)), j = as.integer(factor(long$referenced_works)),
                    x = 1, dimnames = list(levels(factor(long$doi)), NULL))
  association_strength(B)                                    # same normalization as WoS / Scopus
}
adj_oa  <- coupling_oa(corpus_oa)
comm_oa <- louvain_communities(adj_oa, corpus_oa$doi) |> rename(doi = id)

## ---- Etiquetado de comunidades (semiautomático, revisión manual obligatoria) ----
# Bare "invasion" removed: it captures trophoblast invasion (placenta, preeclampsia). Tumour invasion is still
# caught by cancer/tumour/carcinoma/metastasis in the same abstract.
onc_regex   <- regex("cancer|tumou?r|carcinoma|metasta|oncolog|melanoma|neoplas", ignore_case = TRUE)
# Reproductive papers introduce KISS1 as "originally identified as a metastasis suppressor gene": the phrase is
# removed before matching so it does not count as oncology vocabulary (real cancer papers keep their other terms)
onc_boilerplate <- regex("(a |the )?(human )?metastasis[- ]suppress[a-z]*( gene)?", ignore_case = TRUE)
neuro_regex <- regex("GnRH|puberty|hypothalam|LH secretion|gonadotropin|KNDy|neuroendocrin|fertility|ovulat|reproduct", ignore_case = TRUE)

# Neuroendocrine vocabulary ("reproduction", "fertility", "gonadotropin") appears in nearly every kisspeptin
# abstract, cancer papers included, so absolute shares do not separate communities: lift against the
# corpus-wide rate is reported next to them.
# Rule: small clusters -> other; clear oncology enrichment -> oncology; otherwise neuroendocrine (the field's
# default vocabulary) unless that vocabulary is itself depleted. Thresholds live in params (preregistered).
assign_community <- function(n, p_onc, p_neuro, lift_onc, lift_neuro) {
  case_when(
    n < params$community_min_n            ~ "other",
    lift_onc >= params$onc_lift_cut       ~ "oncology",
    lift_neuro >= params$neuro_lift_floor ~ "neuroendocrine",
    TRUE                                  ~ "other")
}

label_clusters <- function(comm, meta) {
  comm |>
    left_join(meta, by = "doi") |>
    mutate(onc   = str_detect(str_remove_all(text_ta, onc_boilerplate), onc_regex),
           neuro = str_detect(text_ta, neuro_regex)) |>
    mutate(base_onc = mean(onc), base_neuro = mean(neuro)) |>            # corpus-wide rates (all labelled papers)
    group_by(cluster) |>
    summarise(n = n(), p_onc = mean(onc), p_neuro = mean(neuro),
              lift_onc = p_onc / first(base_onc), lift_neuro = p_neuro / first(base_neuro), .groups = "drop") |>
    mutate(community_auto = assign_community(n, p_onc, p_neuro, lift_onc, lift_neuro))
}

# Most enriched author keywords per cluster, to support the manual reading of each community
top_keywords <- function(comm, M, k = 10, min_n = 8) {
  kw <- comm |> inner_join(M |> transmute(doi = DI, DE = str_to_lower(coalesce(DE, ""))), by = "doi") |>
    separate_rows(DE, sep = ";") |> mutate(DE = str_trim(DE)) |> filter(DE != "") |> distinct(cluster, doi, DE)
  N <- n_distinct(kw$doi); tot <- count(kw, DE, name = "n_all")
  kw |> count(cluster, DE) |> inner_join(tot, by = "DE") |> filter(n >= min_n) |>
    group_by(cluster) |> mutate(lift = (n / n_distinct(kw$doi[kw$cluster == first(cluster)])) / (n_all / N)) |>
    slice_max(lift, n = k, with_ties = FALSE) |>
    summarise(top_keywords = paste0(DE, " (", round(lift, 1), ")", collapse = "; "), .groups = "drop")
}

meta_main <- M_main |> transmute(doi = DI, year = PY,
                                 text_ta = str_c(coalesce(TI, ""), " ", coalesce(AB, "")))   # NA title would make the label NA
meta_oa   <- corpus_oa |> filter(!is.na(doi)) |> distinct(doi, .keep_all = TRUE) |>        # duplicated DOIs would double-count
  transmute(doi, year = publication_year, text_ta = str_c(coalesce(title, ""), " ", coalesce(abstract, "")))

## ---- Segundo nivel: separar oncología de biología molecular del receptor ----
# The first-level "oncology" cluster also gathers receptor/molecular papers that cite the same early literature
# (Lee 1996, Ohtaki 2001, Kotani 2001). Each oncology cluster is re-partitioned with Louvain on its own
# coupling subgraph and the sub-clusters are relabelled with the same rule (lift against the whole corpus).
refine_oncology <- function(comm, adj, labels) {
  comm$parent_cluster <- comm$cluster
  targets <- labels$cluster[labels$community_auto == "oncology" & labels$n >= 2 * params$community_min_n]
  next_id <- max(comm$cluster)
  for (k in targets) {
    ids <- intersect(comm$doi[comm$cluster == k], rownames(adj))
    sub <- louvain_communities(adj[ids, ids], ids)
    if (n_distinct(sub$cluster) < 2) next
    hit <- match(comm$doi, sub$id)
    comm$cluster[!is.na(hit)] <- next_id + sub$cluster[hit[!is.na(hit)]]   # isolates in the subgraph keep id k
    next_id <- max(comm$cluster)
  }
  comm
}

labels_main <- label_clusters(comm_main, meta_main)
labels_oa   <- label_clusters(comm_oa, meta_oa)
comm_main   <- refine_oncology(comm_main, adj_main, labels_main)
comm_oa     <- refine_oncology(comm_oa, adj_oa, labels_oa)

labels_main <- label_clusters(comm_main, meta_main) |>
  left_join(distinct(comm_main, cluster, parent_cluster), by = "cluster") |>
  left_join(top_keywords(comm_main, M_main), by = "cluster")
labels_oa   <- label_clusters(comm_oa, meta_oa)
write_csv(labels_main, file.path(paths$tables, "community_labels_main.csv"))

## Manual confirmation (main base): data/annotation/community_labels_manual.csv with columns
## cluster, community (oncology / neuroendocrine / other), subfield, reviewer. It overrides community_auto.
manual_file <- file.path(paths$annotation, "community_labels_manual.csv")
if (file.exists(manual_file)) {
  manual <- read_csv(manual_file, col_types = cols(cluster = "i", n = "i", .default = "c"))
  stopifnot(setequal(manual$cluster, labels_main$cluster), all(!is.na(manual$reviewer)))
  # Cluster ids are only meaningful for the partition they were read on: if sizes differ, upstream data changed
  # and Louvain renumbered, so labels must be re-confirmed instead of landing on the wrong groups
  size_check <- inner_join(select(manual, cluster, n_manual = n), select(labels_main, cluster, n), by = "cluster")
  if (any(size_check$n_manual != size_check$n))
    stop("partition changed since manual labelling (cluster sizes differ): re-confirm ", manual_file, call. = FALSE)
  labels_main <- labels_main |> left_join(select(manual, cluster, community_manual = community), by = "cluster")
} else {
  warning("community labels not confirmed by hand yet: using community_auto (", manual_file, ")", call. = FALSE)
  labels_main$community_manual <- labels_main$community_auto
}

comm_main <- comm_main |> left_join(select(labels_main, cluster, community = community_manual), by = "cluster")
comm_oa   <- comm_oa   |> left_join(select(labels_oa,   cluster, community = community_auto), by = "cluster")

## ---- Índice de conectividad entre comunidades por quinquenio ----
connectivity_index <- function(adj, comm, meta) {
  lab    <- c("oncology", "neuroendocrine")
  period <- function(y) cut(y, breaks = seq(1995, 2025, by = 5), labels = FALSE)
  # Composition per period. If coupling ignored community, the expected cross-community share of weight
  # would be 2 p (1 - p), p = oncology share of papers: the raw share moves with p even without any change
  # in behaviour, so obs_exp (observed / expected) is the quantity read for H1.
  comp <- comm |> filter(community %in% lab) |> left_join(select(meta, doi, year), by = "doi") |>
    mutate(period = period(year)) |>
    group_by(period) |> summarise(n_papers = n(), p_onc = mean(community == "oncology"), .groups = "drop")
  ed  <- summary(as(adj, "dgCMatrix")) |> as_tibble() |> filter(i < j)
  ids <- rownames(adj)
  ed |>
    mutate(a = ids[i], b = ids[j]) |>
    left_join(select(comm, doi, ca = community), by = c("a" = "doi")) |>
    left_join(select(comm, doi, cb = community), by = c("b" = "doi")) |>
    left_join(select(meta, doi, year), by = c("a" = "doi")) |>
    filter(ca %in% lab, cb %in% lab) |>
    mutate(period = period(year)) |>
    group_by(period) |>
    summarise(n_edges = n(), cross_share = sum(x[ca != cb]) / sum(x), .groups = "drop") |>  # cross-community weight / total
    left_join(comp, by = "period") |>
    mutate(years    = str_c(1991 + 5 * period, "-", 1995 + 5 * period),
           expected = 2 * p_onc * (1 - p_onc),
           obs_exp  = cross_share / expected) |>
    relocate(period, years)
}
conn_main <- connectivity_index(adj_main, comm_main, meta_main)
write_csv(conn_main, file.path(paths$tables, "connectivity_by_period.csv"))

## ---- Co-citación (base intelectual) y mapa temático ----
cocit <- biblioNetwork(M_main, analysis = "co-citation", network = "references", sep = ";")
saveRDS(cocit, file.path(paths$processed, "cocitation_matrix.rds"))
tm <- thematicMap(M_main, field = "DE", n = 250, minfreq = 5)           # author keywords
saveRDS(tm, file.path(paths$processed, "thematic_map.rds"))
# Burst detection on yearly keyword counts: bursts::kleinberg() per keyword (verify package availability)

## ---- Sensibilidad entre bases: ARI sobre DOIs comunes ----
ari_tbl <- inner_join(select(comm_main, doi, c_main = cluster),
                      select(comm_oa,   doi, c_oa   = cluster), by = "doi")
ari <- tibble(
  n_common  = nrow(ari_tbl),
  ari       = mclust::adjustedRandIndex(ari_tbl$c_main, ari_tbl$c_oa),
  ari_label = mclust::adjustedRandIndex(                 # agreement on the substantive labels (H1)
    comm_main$community[match(ari_tbl$doi, comm_main$doi)],
    comm_oa$community[match(ari_tbl$doi, comm_oa$doi)])
)
print(ari)
write_csv(ari, file.path(paths$supp, "S2_ari_between_bases.csv"))

save_versioned(comm_main, "communities_main")
save_versioned(comm_oa,   "communities_openalex")
