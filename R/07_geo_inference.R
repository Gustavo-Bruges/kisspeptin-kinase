## =============================================================================
## 07_geo_inference.R · Rama de datos: limma, decoupleR, CARNIVAL + controles
## REGLA DE INDEPENDENCIA: este script no lee salidas de 03–06
## =============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(decoupleR); library(httr2)
  library(CARNIVAL)
})

## ---- Series elegibles y contrastes (definidos sobre el pData real; decisions.md 2026-10-18) ----
# GEO and ArrayExpress were searched (2026-10-18): no third series perturbs KISS1R with n >= 3 per group,
# so the backup criterion of PROYECTO.md applies (data module presented as exploratory).
series <- tribble(
  ~gse,       ~organism, ~sign, ~note,
  "GSE25458", "human",   +1,    "Ishikawa endometrial cancer cells; 2x2 siRNA (ctrl / GPR54) x peptide (ctrl / metastin-10), n = 3",
  "GSE28383", "mouse",   -1,    "male hypothalamus; GPR54-/- (4) + KISS1-/- (5) vs WT 3 wk (3) and WT 6 wk (5); 3 Affy batches"
)

## Expression, probe symbols and sample annotation per series
first_symbol <- function(x, sep) str_trim(str_split_i(x, sep, 1))
load_series <- list(
  GSE25458 = function() {
    e  <- getGEO("GSE25458", GSEMatrix = TRUE, destdir = paths$raw_geo)[[1]]      # GPL570, RMA log2 in series matrix
    pd <- pData(e)
    list(ex = exprs(e), symbol = first_symbol(fData(e)[["Gene Symbol"]], "///"),
         pd = tibble(sample  = rownames(pd),
                     sirna   = if_else(str_detect(pd[["sirna transfection:ch1"]], "GPR54"), "siGPR54", "siCtrl"),
                     peptide = if_else(str_detect(pd[["peptide treatment:ch1"]], "metastin"), "kp10", "ctrl")) |>
           mutate(group = str_c(sirna, "_", peptide)))
  },
  GSE28383 = function() {
    # Series matrix has 0 rows (only CEL files deposited): expression comes from R/07a_process_cel_GSE28383.R
    obj <- readRDS(file.path(paths$raw_geo, "GSE28383_rma_core.rds"))
    pd  <- obj$pd
    list(ex = obj$ex, symbol = obj$symbol,
         pd = tibble(sample = pd$sample,
                     group  = case_when(str_detect(pd$genotype, "-/-") ~ "KO",          # GPR54-/- and KISS1-/- pooled
                                        str_detect(pd$age, "3 weeks")   ~ "WT3",
                                        str_detect(pd$age, "6 weeks")   ~ "WT6"),
                     batch  = factor(str_extract(pd$source_name, "batch [0-9]+"))))     # Affy processing batch: covariate
  }
)

## Fits per series: each fit has its own sample subset, design and contrasts. Contrast roles:
##   "network"        CARNIVAL + stability; the series signal is the intersection over its network contrasts
##   "receptor_check" single CARNIVAL run; kinases that persist flag receptor-independent responses (GSE25458)
##   "sensitivity"    CARNIVAL + stability, reported separately, never part of the series signal
## GSE28383: processing batch 2 (4 arrays) dominates PC1 (60 % of variance) with a compressed distribution
## (IQR 2.1 vs 3.0): main fit keeps all arrays with batch covariate + limma::arrayWeights; sensitivity drops
## batch 2 (WT3 would keep n = 2, so only KO vs WT6 survives). decisions.md 2026-10-18.
contrast_defs <- list(
  GSE25458 = list(
    main = list(subset = NULL, design = ~ 0 + group, array_weights = FALSE,
                contrasts = c(kp10_siCtrl  = "groupsiCtrl_kp10 - groupsiCtrl_ctrl",
                              kp10_siGPR54 = "groupsiGPR54_kp10 - groupsiGPR54_ctrl"),
                role = c(kp10_siCtrl = "network", kp10_siGPR54 = "receptor_check"))),
  GSE28383 = list(
    main = list(subset = NULL, design = ~ 0 + group + batch, array_weights = TRUE,
                contrasts = c(KO_vs_WT6 = "groupKO - groupWT6",       # loss of KISS1R + hypogonadism
                              KO_vs_WT3 = "groupKO - groupWT3"),      # vs prepubertal WT: partly hormone-matched
                role = c(KO_vs_WT6 = "network", KO_vs_WT3 = "network")),
    no_batch2 = list(subset = \(pd) pd$batch != "batch 2" & pd$group %in% c("KO", "WT6"),
                     design = ~ 0 + group + batch, array_weights = FALSE,
                     contrasts = c(KO_vs_WT6_noB2 = "groupKO - groupWT6"),
                     role = c(KO_vs_WT6_noB2 = "sensitivity")))
)

## ---- Expresión diferencial por ajuste ----
run_de <- function(dat, fit_def) {
  keep <- if (is.null(fit_def$subset)) rep(TRUE, nrow(dat$pd)) else fit_def$subset(dat$pd)
  pd <- droplevels(dat$pd[keep, ]); ex <- dat$ex[, keep, drop = FALSE]
  stopifnot(identical(colnames(ex), pd$sample), min(table(pd$group)) >= 3)   # eligibility: n >= 3 per group
  if (max(ex, na.rm = TRUE) > 100) ex <- log2(ex + 1)        # series-matrix values sometimes unlogged
  X  <- model.matrix(fit_def$design, data = pd)
  colnames(X) <- make.names(colnames(X))                     # "batchbatch 2" is not a valid contrast level
  aw <- if (isTRUE(fit_def$array_weights)) arrayWeights(ex, X) else NULL      # down-weights noisier arrays
  fit <- lmFit(ex, X, weights = aw)
  cm  <- makeContrasts(contrasts = fit_def$contrasts, levels = X); colnames(cm) <- names(fit_def$contrasts)
  fit2 <- contrasts.fit(fit, cm) |> eBayes()
  tts <- imap(fit_def$contrasts, \(x, k) {
    topTable(fit2, coef = k, number = Inf, sort.by = "none") |> rownames_to_column("probe") |>
      mutate(symbol = dat$symbol[match(probe, rownames(dat$ex))]) |>
      filter(!is.na(symbol), symbol != "", symbol != "---") |>
      group_by(symbol) |> slice_max(AveExpr, n = 1, with_ties = FALSE) |> ungroup()  # one probe per gene: highest mean expression, blind to the contrast
  })
  attr(tts, "array_weights") <- if (is.null(aw)) NULL else setNames(aw, colnames(ex))
  tts
}

## ---- OmniPath por la API web (OmnipathR 4.0.0 falla en el cruce por ncbi_tax_id; decisions.md 2026-10-18) ----
# Same source as OmnipathR, cached as dated TSV in data/raw/curated. Logical fields come as "True"/"False".
fetch_omnipath <- function(name, query) {
  f <- sort(list.files(paths$curated, sprintf("^omnipath_%s_.*[.]tsv$", name), full.names = TRUE), decreasing = TRUE)[1]
  if (is.na(f)) {
    f <- file.path(paths$curated, sprintf("omnipath_%s_%s.tsv", name, Sys.Date()))
    request("https://omnipathdb.org/interactions") |>
      req_url_query(!!!query, genesymbols = "yes", format = "tsv", license = "academic") |>
      req_retry(max_tries = 3) |> req_perform(path = f)
  }
  read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c")) |>
    mutate(across(c(is_directed, is_stimulation, is_inhibition, consensus_direction,
                    consensus_stimulation, consensus_inhibition), \(x) x == "True"))
}

## CollecTRI regulons in decoupleR format (source, target, mor); sign as decoupleR::get_collectri():
## inhibition-only -> -1, otherwise +1. Complexes keep their OmniPath genesymbol (e.g. "FOS_JUN").
get_collectri_web <- function(organism) {
  fetch_omnipath(str_c("collectri_", organism),
                 list(datasets = "collectri", organisms = c(human = 9606, mouse = 10090)[[organism]])) |>
    transmute(source = source_genesymbol, target = target_genesymbol,
              mor = if_else(consensus_inhibition & !consensus_stimulation, -1, 1)) |>
    distinct(source, target, .keep_all = TRUE)
}

## PROGENy from the Bioconductor package model (decoupleR::get_progeny downloads it through the broken OmnipathR):
## same selection as get_progeny(top = 500), the 500 genes with lowest p-value per pathway
get_progeny_local <- function(organism, top = 500) {
  m <- if (organism == "mouse") progeny::model_mouse_full else progeny::model_human_full
  as_tibble(m) |> group_by(pathway) |> slice_min(p.value, n = top, with_ties = FALSE) |> ungroup() |>
    transmute(source = as.character(pathway), target = as.character(gene), weight, p_value = p.value)   # factors in the package
}

## ---- Actividad de TF y vías ----
run_activity <- function(tt, organism) {
  mat <- matrix(tt$t, dimnames = list(tt$symbol, "t"))      # t-statistic keeps direction and precision
  list(
    tf = run_ulm(mat, get_collectri_web(organism),
                 .source = "source", .target = "target", .mor = "mor", minsize = 5),
    pw = run_mlm(mat, get_progeny_local(organism, top = 500),
                 .source = "source", .target = "target", .mor = "weight", minsize = 5)
  )
}

## ---- Red previa OmniPath con signo ----
pkn_full <- fetch_omnipath("interactions", list(fields = "references,curation_effort")) |>
  filter(is_directed, consensus_stimulation != consensus_inhibition) |>   # unambiguously signed only
  transmute(source = source_genesymbol,
            interaction = if_else(consensus_stimulation, 1L, -1L),
            target = target_genesymbol,
            references)

## Control de circularidad: quitar aristas sostenidas solo por artículos del corpus de kisspeptina
# corpus_pmids.txt is exported from 02 as a plain list, so this script never reads analysis outputs of 03–06
corpus_pmids <- readLines(file.path(paths$processed, "corpus_pmids.txt"))
edge_only_kiss <- function(refs) {
  pm <- str_extract_all(refs, "\\d{6,9}")[[1]]
  length(pm) > 0 && all(pm %in% corpus_pmids)
}
pkn_nocirc <- pkn_full |> filter(!map_lgl(references, edge_only_kiss))

## Fake receptor (negative control), fixed by rule (decisions.md 2026-10-18): GPCR whose PKN edges reach Gs/Gi
## but no Gq-family protein (KISS1R's only edges are GNAQ and GNA14, so a Gq receptor would reproduce its network),
## AveExpr below the 25th percentile in both series, out-degree equal to KISS1R's; ties -> lowest expression.
## HCRTR1, the earlier proposal, was dropped: its 8 G-protein edges include both KISS1R edges.
fake_receptor <- "HRH4"
gq_family <- c("GNAQ", "GNA11", "GNA14", "GNA15")
# Perturbed nodes need outgoing edges, otherwise CARNIVAL cannot propagate from them
stopifnot(all(c("KISS1R", fake_receptor) %in% pkn_nocirc$source),
          !any(pkn_nocirc$target[pkn_nocirc$source == fake_receptor] %in% gq_family))

# Node degree in the PKN actually used, for the hub-controlled transfer test in 08
pkn_degree <- tibble(symbol = c(pkn_nocirc$source, pkn_nocirc$target)) |> count(symbol, name = "degree")
saveRDS(pkn_degree, file.path(paths$processed, "pkn_nocirc_degree.rds"))

## ---- CARNIVAL ----
to_human_tf <- function(tf_scores, organism) {
  if (organism == "human") return(tf_scores)
  orth <- babelgene::orthologs(genes = tf_scores$source, species = organism, human = FALSE) |>
    select(symbol, human_symbol) |> distinct()
  tf_scores |> inner_join(orth, by = c("source" = "symbol")) |> mutate(source = human_symbol)
}

## Solver (params$carnival_solver). lpSolve builds a dense constraint matrix (asked for 700 GB on the full PKN).
##   "gurobi": academic licence, gurobi_cl on PATH or GUROBI_CL env var; needed for depth >= 6
##   "cbc":    CBC 2.10.12 (COIN-OR, open source) outside the R project, or CBC_PATH env var; optimum only up to depth 4
##   "highs":  HiGHS 1.15.1 (MIT, open source) outside the R project, or HIGHS_PATH env var; plugged into CARNIVAL's
##             CBC slot (CARNIVAL writes the LP file, HiGHS solves it, the solution is returned in CBC's layout)
solver_path <- switch(params$carnival_solver,
  gurobi = Sys.getenv("GUROBI_CL", unname(Sys.which("gurobi_cl"))),
  highs  = Sys.getenv("HIGHS_PATH", file.path(dirname(here::here()), "tools", "highs-1.15.1", "bin", "highs.exe")),
  cbc    = Sys.getenv("CBC_PATH", file.path(dirname(here::here()), "tools", "cbc-2.10.12", "bin", "cbc.exe")))

## HiGHS raw solution file: "Model status" / <status>, then "# Columns n" followed by "<name> <value>" lines
solve_highs <- function(carnivalOptions) {
  base <- str_remove(carnivalOptions$filenames$resultFile, "[.]txt$")
  opt  <- str_c(base, "_highs.opt"); sol <- str_c(base, "_highs.sol")
  writeLines(c(sprintf("mip_rel_gap = %g", carnivalOptions$poolrelGap),
               sprintf("random_seed = %d", params$seed %% .Machine$integer.max)), opt)
  log <- system2(carnivalOptions$solverPath, c("--model_file", carnivalOptions$filenames$lpFilename,
                 "--options_file", opt, "--time_limit", carnivalOptions$timelimit,
                 "--threads", carnivalOptions$threads, "--solution_file", sol), stdout = TRUE, stderr = TRUE)
  txt <- readLines(sol)
  status <- txt[which(txt == "Model status")[1] + 1]
  i <- grep("^# Columns ", txt)[1]; n <- as.integer(str_extract(txt[i], "[0-9]+$"))
  cols <- str_split_fixed(txt[(i + 1):(i + n)], " ", 2)
  solver_status$last <- list(optimal = identical(status, "Optimal"), result = status,
                             log_tail = tail(log, 15))
  data.frame(name = cols[, 1], solution = as.character(round(as.numeric(cols[, 2]))))   # binaries/integers: drop 1e-9 noise
}
stopifnot(nzchar(solver_path), file.exists(solver_path))

## CARNIVAL 2.22 reads the solution file whatever the solver status: a CBC run stopped on time limit without any
## feasible point (objective 1e+50) still yielded a "network". Same calls as CARNIVAL:::solveWithCbc / solveWithGurobi,
## but the log is captured and the status kept in solver_status$last, checked by run_carnival().
solver_status <- new.env()
assignInNamespace("solveWithCbc", function(carnivalOptions) {
  if (params$carnival_solver == "highs") return(solve_highs(carnivalOptions))
  resultFile <- str_replace(carnivalOptions$filenames$resultFile, ".txt", ".csv")
  log <- system2(carnivalOptions$solverPath, c("-import", carnivalOptions$filenames$lpFilename,
                 "-seconds", carnivalOptions$timelimit, "-ratioGap", carnivalOptions$poolrelGap,
                 "-threads", carnivalOptions$threads, "-solve", "-printi", "csv", "-solution", resultFile),
                 stdout = TRUE, stderr = TRUE)
  res <- str_subset(log, "^Result - ")
  solver_status$last <- list(optimal = any(str_detect(res, "Optimal solution found")),
                             result  = if (length(res)) str_remove(res[1], "^Result - ") else "unknown")
  read.csv2(resultFile, sep = ",")
}, ns = "CARNIVAL")
assignInNamespace("solveWithGurobi", function(carnivalOptions) {
  resultFile <- str_sub(carnivalOptions$filenames$resultFile, end = -5)
  log <- system2(carnivalOptions$solverPath, c(
    paste0("MIPGap=", carnivalOptions$mipGap), paste0("TimeLimit=", carnivalOptions$timelimit),
    paste0("PoolGap=", carnivalOptions$poolrelGap), paste0("SolutionLimit=", carnivalOptions$limitPop),
    paste0("PoolSolutions=", carnivalOptions$poolCap), paste0("Threads=", carnivalOptions$threads),
    paste0("SolFiles=", resultFile), "PoolSearchMode=2", "NumericFocus=2",
    paste0("ResultFile=", resultFile, ".sol"), carnivalOptions$filenames$lpFilename), stdout = TRUE, stderr = TRUE)
  res <- str_subset(log, "^(Optimal solution found|Time limit reached|Model is infeasible|Solution count)")
  solver_status$last <- list(optimal = any(str_detect(log, "^Optimal solution found")),
                             result  = if (length(res)) str_c(res, collapse = "; ") else "unknown")
  resultFile
}, ns = "CARNIVAL")

carnival_opts <- function() {
  o <- defaultCbcSolveCarnivalOptions(solverPath = solver_path)
  if (params$carnival_solver == "gurobi")
    o <- modifyList(o, list(solver = "gurobi", mipGap = 1e-4,
                            limitPop = 2e9,          # gurobi SolutionLimit: stop after N feasible points; effectively off
                            poolCap = 100))          # equally optimal networks kept -> weightedSIF Weight = % of them
  wd <- file.path(tempdir(), "carnival"); dir.create(wd, showWarnings = FALSE)
  o$outputFolder <- wd; o$workdir <- wd; o$keepLPFiles <- FALSE   # defaults write LP files into the project root
  o$timelimit <- params$carnival_timelimit; o$threads <- params$carnival_threads
  o
}

## Path pruning: keep node v only if it lies on a path perturbed node -> v -> measured TF of at most
## params$pkn_max_depth steps, i.e. d(from, v) + min_m d(v, m) <= depth. Nodes on no such path cannot enter
## the solution; they only enlarge the ILP (11 k edges without the depth cap: CBC found no feasible point in 600 s).
prune_pkn <- function(pkn, from, meas_nodes, depth = params$pkn_max_depth) {
  g <- igraph::graph_from_data_frame(distinct(pkn, source, target), directed = TRUE)
  meas_nodes <- intersect(meas_nodes, igraph::V(g)$name)
  if (!length(meas_nodes)) return(pkn[0, ])
  d_from <- igraph::distances(g, v = from, mode = "out")[1, ]
  d_to   <- apply(igraph::distances(g, to = meas_nodes, mode = "out"), 1, min)
  keep   <- names(which(d_from + d_to <= depth))
  filter(pkn, source %in% keep, target %in% keep)
}

run_carnival <- function(tf, sign, pkn, perturb_node = "KISS1R") {
  meas <- tf |> filter(p_value < 0.05) |> (\(d) setNames(d$score, d$source))()
  pkn  <- prune_pkn(pkn, perturb_node, names(meas))
  meas <- meas[names(meas) %in% c(pkn$source, pkn$target)]   # TFs with no path from the perturbed node drop out
  if (!length(meas)) return(tibble(Node1 = character(), Sign = numeric(), Node2 = character(), Weight = numeric()))
  solver_status$last <- NULL                                 # never inherit the status of the previous run
  res <- runVanillaCarnival(
    perturbations = setNames(sign, perturb_node),
    measurements  = meas,
    priorKnowledgeNetwork = as.data.frame(select(pkn, source, interaction, target)),
    carnivalOptions = carnival_opts()
  )
  # weightedSIF lists EVERY PKN edge with Weight = % of solutions containing it (0 = absent). Without any
  # solution (solver stopped on time limit) all weights are 0/NA: fail loudly instead of returning the PKN.
  # Only proven optima are kept: a time-limited incumbent depends on run time and threads, not on the data
  if (!isTRUE(solver_status$last$optimal))
    stop(params$carnival_solver, ": ", solver_status$last$result, " (", perturb_node, ")")
  w <- as_tibble(res$weightedSIF)
  if (!any(coalesce(w$Weight, 0) > 0)) stop("CARNIVAL returned no solution (", perturb_node, ")")
  out <- filter(w, Weight > 0)
  attr(out, "nodes") <- filter(as_tibble(res$nodesAttributes), AvgAct != 0)   # node activity sign (H2 sign check)
  attr(out, "n_meas") <- length(meas); attr(out, "n_pkn_edges") <- nrow(pkn)
  out
}

## Same seeds for every perturbed node: the fake receptor sees exactly the same TF subsamples as KISS1R (paired)
stability_nodes <- function(tf, sign, pkn, perturb_node = "KISS1R") {
  runs <- map(seq_len(params$n_stability), \(i) {
    set.seed(params$seed + i)
    sub <- slice_sample(tf, prop = params$stability_frac)
    sif <- tryCatch(run_carnival(sub, sign, pkn, perturb_node), error = function(e) NULL)
    if (i %% 10 == 0) message(sprintf("[%s] stability %d/%d", format(Sys.time(), "%H:%M"), i, params$n_stability))
    if (is.null(sif)) NULL else unique(c(sif$Node1, sif$Node2))
  })
  ok <- !map_lgl(runs, is.null)            # failed runs (no solution in time) leave the denominator, never count as "absent"
  if (mean(ok) < 0.9) warning(sprintf("stability: only %d of %d CARNIVAL runs solved", sum(ok), length(ok)))
  out <- tibble(node = unlist(runs[ok])) |> count(node) |>
    mutate(freq = n / sum(ok)) |>
    filter(freq >= params$stability_keep)
  attr(out, "n_solved") <- sum(ok)
  out
}

## ---- Ejecución por serie, ajuste y contraste ----
run_contrast <- function(tt, organism, sign, role) {
  act <- run_activity(tt, organism)
  tfh <- to_human_tf(act$tf, organism)
  if (role == "receptor_check")
    return(list(role = role, de = tt, activity = act, carnival_nocirc = run_carnival(tfh, sign, pkn_nocirc)))
  list(role = role, de = tt, activity = act,
       carnival        = run_carnival(tfh, sign, pkn_full),
       carnival_nocirc = run_carnival(tfh, sign, pkn_nocirc),
       carnival_fake   = run_carnival(tfh, sign, pkn_nocirc, perturb_node = fake_receptor),  # same PKN as the real run
       stable_nodes    = stability_nodes(tfh, sign, pkn_nocirc))                           # "stable and non-circular"
}

## Each contrast is cached as soon as it finishes (a stability block takes hours at depth 5). The file name
## carries solver and depth, so changing either never reuses stale networks. KISS_GSE=<gse> runs one series only
## (series in parallel, one R process each); a final run without it finds every cache and builds data_kin.
cache_dir <- file.path(paths$processed, "carnival_cache"); dir.create(cache_dir, showWarnings = FALSE)
only_gse  <- Sys.getenv("KISS_GSE")
run_series <- if (nzchar(only_gse)) filter(series, gse == only_gse) else series
stopifnot(nrow(run_series) > 0)

results <- pmap(run_series, \(gse, organism, sign, note) {
  dat  <- load_series[[gse]]()
  fits <- contrast_defs[[gse]]
  per_contrast <- list(); weights <- list()
  for (f in names(fits)) {
    des <- run_de(dat, fits[[f]])
    weights[[f]] <- attr(des, "array_weights")
    for (k in names(des)) {
      cf <- file.path(cache_dir, sprintf("%s_%s_%s_d%d_n%d.rds", gse, k, params$carnival_solver,
                                         params$pkn_max_depth, params$n_stability))
      if (!file.exists(cf)) {
        message(sprintf("[%s] %s %s start", format(Sys.time(), "%H:%M"), gse, k))
        saveRDS(run_contrast(des[[k]], organism, sign, fits[[f]]$role[[k]]), cf)
      }
      per_contrast[[k]] <- readRDS(cf)
      # Fake receptor under the same stability protocol as KISS1R (a single fake run vs 100 real subsamples would
      # make the control look more permissive than it is). Own cache: the real stability block is never recomputed.
      if (per_contrast[[k]]$role != "receptor_check") {
        ff <- str_replace(cf, "[.]rds$", sprintf("_fake_%s.rds", fake_receptor))
        if (!file.exists(ff)) {
          message(sprintf("[%s] %s %s fake stability start", format(Sys.time(), "%H:%M"), gse, k))
          tfh <- to_human_tf(per_contrast[[k]]$activity$tf, organism)
          saveRDS(stability_nodes(tfh, sign, pkn_nocirc, perturb_node = fake_receptor), ff)
        }
        per_contrast[[k]]$stable_nodes_fake <- readRDS(ff)
      }
    }
  }
  roles <- map_chr(per_contrast, "role")
  series_signal <- \(slot) reduce(map(per_contrast[roles == "network"], \(x) x[[slot]]$node), intersect)
  list(gse = gse, note = note, contrasts = per_contrast, array_weights = weights,
       stable_nodes      = series_signal("stable_nodes"),
       stable_nodes_fake = series_signal("stable_nodes_fake"))     # same rule as the real signal
})
names(results) <- run_series$gse
if (nzchar(only_gse)) { message("KISS_GSE=", only_gse, ": cached; rerun without KISS_GSE to build data_kin"); quit(save = "no") }
saveRDS(results, file.path(paths$processed, "geo_inference_results.rds"))

## ---- Soporte en datos por quinasa ----
kinome <- read_csv(file.path(paths$curated, "kinome_reference.csv"), show_col_types = FALSE)$symbol
sif_nodes <- function(sif) unique(c(sif$Node1, sif$Node2))
by_role   <- function(r, role) keep(r$contrasts, \(x) x$role == role)
data_kin <- imap_dfr(results, \(r, gse) {
  net <- by_role(r, "network"); check <- by_role(r, "receptor_check"); sens <- by_role(r, "sensitivity")
  tibble(symbol = intersect(r$stable_nodes, kinome), gse = gse) |>
    mutate(in_fake      = symbol %in% r$stable_nodes_fake,                              # stable under HRH4 too
           in_fake_single = symbol %in% unlist(map(net, \(x) sif_nodes(x$carnival_fake))),  # single fake run (descriptive)
           receptor_ind = if (length(check)) symbol %in% unlist(map(check, \(x) sif_nodes(x$carnival_nocirc))) else NA,
           in_sens      = if (length(sens)) symbol %in% unlist(map(sens, \(x) x$stable_nodes$node)) else NA)
}) |>
  group_by(symbol) |>
  summarise(data_freq = n_distinct(gse) / nrow(series),
            in_fake_any = any(in_fake), in_fake_single_any = any(in_fake_single),
            receptor_independent = if (all(is.na(receptor_ind))) NA else any(receptor_ind, na.rm = TRUE),
            survives_sensitivity = if (all(is.na(in_sens))) NA else all(in_sens, na.rm = TRUE),  # NA: no sensitivity run
            .groups = "drop")
save_versioned(data_kin, "data_kin")
write_csv(data_kin, file.path(paths$tables, "data_kin.csv"))
