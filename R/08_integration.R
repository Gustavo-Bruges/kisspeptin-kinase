## =============================================================================
## 08_integration.R · Cruce literatura × datos (H2, H3, análisis temporal)
## Umbrales tomados de params (preregistrados)
## =============================================================================
source(here::here("R", "00_setup.R"))

lit_kin  <- load_versioned("lit_kin")
data_kin <- load_versioned("data_kin")
results  <- readRDS(file.path(paths$processed, "geo_inference_results.rds"))

# Canonical KISS1R kinases named in H2 (Gq -> PKC -> ERK1/2; p38; PI3K/AKT -| GSK3B)
canonical <- c("MAPK1", "MAPK3", "MAPK14", "AKT1", "PIK3CA", "GSK3B")

## ---- Cuadrantes ----
# "High literature support" needs >= params$lit_min_docs corpus papers: with 1-2 papers and a small PubMed
# background, lit_score is noise. The p75 cut is also computed only over eligible kinases.
lit_eligible <- lit_kin$n_docs >= params$lit_min_docs
lit_cut <- quantile(lit_kin$lit_score[lit_eligible], params$lit_high_quantile, na.rm = TRUE)

integrated <- full_join(lit_kin, data_kin, by = "symbol") |>
  mutate(
    lit_score = coalesce(lit_score, min(lit_kin$lit_score, na.rm = TRUE)),  # kinase absent from literature = lowest support
    data_freq = coalesce(data_freq, 0),
    lit_high  = lit_score > lit_cut & coalesce(n_docs, 0) >= params$lit_min_docs,
    data_high = data_freq >= params$data_high_frac,
    canonical = symbol %in% canonical,
    # specific = stable with KISS1R but not under the fake receptor (same stability protocol) and gone when
    # KISS1R is silenced; NA when the series has no siGPR54 arm
    specific  = data_high & !coalesce(in_fake_any, FALSE) & !coalesce(receptor_independent, FALSE),
    quadrant  = case_when(
      lit_high &  data_high ~ "confirmed",
      !lit_high & data_high ~ "novel_hypothesis",
      lit_high & !data_high ~ "not_reproduced",
      TRUE                  ~ "no_evidence")
  )
write_csv(integrated, file.path(paths$tables, "kinase_quadrants.csv"))

## Table 2: prioritized hypotheses. Canonical kinases can land in "novel_hypothesis" only because the
## kisspeptin corpus under-cites them (e.g. MAPK14: lit_score below the p75 cut), not because they are new
## to KISS1R signalling.
## Canonical kinases leave Table 2 and are reported apart as "canonical, under-cited in the kisspeptin corpus".
## Order: specific (survives fake receptor and siGPR54) first, then support across series, then the least cited.
prioritize <- function(integrated) {
  integrated |>
    filter(quadrant == "novel_hypothesis", !canonical) |>
    arrange(desc(coalesce(specific, FALSE)), desc(data_freq), lit_score)
}
write_csv(prioritize(integrated), file.path(paths$tables, "table2_prioritized_hypotheses.csv"))
write_csv(filter(integrated, quadrant == "novel_hypothesis", canonical) |> arrange(desc(data_freq)),
          file.path(paths$tables, "canonical_undercited.csv"))

## ---- H2: canonical kinases, KISS1R vs fake receptor under the same stability protocol ----
h2 <- imap_dfr(results, \(r, gse) {
  tibble(gse = gse,
         canonical_real = sum(canonical %in% r$stable_nodes),
         canonical_fake = sum(canonical %in% r$stable_nodes_fake),
         real_only      = str_c(sort(setdiff(intersect(canonical, r$stable_nodes), r$stable_nodes_fake)), collapse = ","),
         fake_only      = str_c(sort(setdiff(intersect(canonical, r$stable_nodes_fake), r$stable_nodes)), collapse = ","),
         n_canonical    = length(canonical))
})
print(h2)
write_csv(h2, file.path(paths$tables, "h2_canonical_recovery.csv"))

## H2 sign check: expected activity under the agonist (sign +1); a knockout (sign -1) flips it.
## GSK3B is inhibited by AKT phosphorylation, so it is expected down under the agonist.
expected_sign <- c(MAPK1 = 1, MAPK3 = 1, MAPK14 = 1, AKT1 = 1, PIK3CA = 1, GSK3B = -1)
series_sign   <- c(GSE25458 = 1, GSE28383 = -1)
h2_sign <- imap_dfr(results, \(r, gse) {
  s <- series_sign[[gse]]                          # outside transmute: there `gse` is already the new column
  imap_dfr(keep(r$contrasts, \(x) x$role == "network"), \(x, k) {
    attr(x$carnival_nocirc, "nodes") |>
      filter(Node %in% canonical) |>
      transmute(gse = gse, contrast = k, symbol = Node, avg_activity = AvgAct,
                expected = expected_sign[symbol] * s,
                sign_ok  = sign(avg_activity) == expected)
  })
})
print(h2_sign)
write_csv(h2_sign, file.path(paths$tables, "h2_sign_check.csv"))

## ---- H3: prueba de transferencia (quinasas oncológicas en la red hipotalámica) ----
hyp_nodes  <- results$GSE28383$stable_nodes          # stable in both KO vs WT6 and KO vs WT3
pkn_degree <- readRDS(file.path(paths$processed, "pkn_nocirc_degree.rds"))

# Only kinases present in the PKN can ever enter a CARNIVAL network; degree bins control for hub bias
# (oncology kinases have more curated edges and would appear in any network)
kin_lab <- lit_kin |>
  filter(community %in% c("oncology", "neuroendocrine")) |>
  inner_join(pkn_degree, by = "symbol") |>
  mutate(deg_bin = ntile(degree, 4)) |>
  select(symbol, community, degree, deg_bin)

stat_transfer <- function(lab) sum(kin_lab$symbol[lab == "oncology"] %in% hyp_nodes)

# Shuffle community labels while preserving how many oncology labels each degree bin holds
permute_labels <- function(lab, strata) {
  out <- lab
  for (s in unique(strata)) {
    i <- which(strata == s)
    out[i] <- lab[i][sample.int(length(i))]                 # sample.int: safe for strata of size 1
  }
  out
}

obs <- stat_transfer(kin_lab$community)
set.seed(params$seed)
null <- replicate(params$n_perm, stat_transfer(permute_labels(kin_lab$community, kin_lab$deg_bin)))
transfer <- tibble(n_hyp_kinases = sum(kin_lab$symbol %in% hyp_nodes),
                   n_oncology = sum(kin_lab$community == "oncology"), n_labelled = nrow(kin_lab),
                   observed = obs, null_mean = mean(null),
                   p_value = (sum(null >= obs) + 1) / (params$n_perm + 1))
print(transfer)
write_csv(transfer, file.path(paths$tables, "h3_transfer_test.csv"))
saveRDS(null, file.path(paths$processed, "transfer_null.rds"))

## ---- Análisis temporal ----
temporal <- integrated |>
  filter(quadrant %in% c("confirmed", "novel_hypothesis"), is.finite(first_year)) |>
  group_by(quadrant) |>
  summarise(n = n(), median_first_year = median(first_year), .groups = "drop")
write_csv(temporal, file.path(paths$tables, "temporal_first_mention.csv"))

## ---- Figura 6: dispersión lit_score × soporte en datos ----
library(ggplot2)
p6 <- ggplot(integrated, aes(lit_score, data_freq, colour = quadrant, shape = canonical)) +
  geom_hline(yintercept = params$data_high_frac, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = lit_cut, linetype = 2, colour = "grey50") +
  geom_jitter(width = 0, height = 0.02, alpha = 0.8, size = 2) +
  geom_text(data = filter(integrated, data_high), aes(label = symbol), size = 3, vjust = -0.8,
            show.legend = FALSE) +
  labs(x = "Soporte en literatura (lit_score, log2 enriquecimiento)",
       y = "Fracción de series con la quinasa en la red causal estable",
       colour = NULL, shape = "Canónica (H2)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(paths$figures, "fig6_quadrants.pdf"), p6, width = 6.5, height = 5)
