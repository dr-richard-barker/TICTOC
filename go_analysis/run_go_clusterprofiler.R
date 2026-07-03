#!/usr/bin/env Rscript
# =============================================================================
# TICTOC — GO / KEGG enrichment of cotton DEGs in Arabidopsis space
# =============================================================================
# Consumes:
#   - DEG tables         TICTOC_3_factor_model/Diff_genes_heatmap_*.csv  (Gohir IDs)
#   - ortholog crosswalk crosswalk/gohir_to_arabidopsis.tsv             (Gohir -> AT)
#   - count matrix       TICTOC_run1_filteredCounts_v3.csv              (for the universe)
#
# Why map to Arabidopsis: cotton (Gohir) has no curated GO/KEGG annotation, but
# Arabidopsis (org.At.tair.db, KEGG "ath") is richly annotated. We translate each
# DEG's Gohir ID to its best-hit AT locus, then enrich in AT space. This is the
# same bridge PhysioSpace will use (roadmap 4.4 / 4.5).
#
# IMPORTANT caveats (read before trusting results):
#   * Best-BLASTP-hit orthology is approximate; many Gohir -> one AT (A/D homoeologs
#     + paralogs) so gene counts COLLAPSE on mapping. Enrichment is over UNIQUE AT loci.
#   * Always enrich against the correct UNIVERSE = AT loci reachable from the EXPRESSED
#     cotton genes (built here from the count matrix), not the whole genome. This script
#     does that; do not drop the `universe=` argument.
#
# Usage:
#   Rscript run_go_clusterprofiler.R [--deg-glob "../TICTOC_3_factor_model/Diff_genes_heatmap_*.csv"]
#                                    [--crosswalk ../crosswalk/gohir_to_arabidopsis.tsv]
#                                    [--counts ../TICTOC_run1_filteredCounts_v3.csv]
#                                    [--out results] [--padj 0.05] [--lfc 1] [--min-pid 30]
#                                    [--use-regulation TRUE] [--kegg TRUE]
#
# Install (once):
#   install.packages("BiocManager")
#   BiocManager::install(c("clusterProfiler","org.At.tair.db","enrichplot"))
#   install.packages(c("ggplot2"))            # optplot rendering
# =============================================================================

suppressWarnings(suppressMessages({
  ok <- requireNamespace("clusterProfiler", quietly = TRUE) &&
        requireNamespace("org.At.tair.db",  quietly = TRUE)
}))
if (!ok) stop("Missing packages. Run:\n  install.packages('BiocManager')\n  ",
              "BiocManager::install(c('clusterProfiler','org.At.tair.db','enrichplot'))")

# ---- tiny CLI parser (no external deps) -------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[[i + 1]]
}
HERE       <- tryCatch(dirname(sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), error = function(e) ".")
if (is.na(HERE) || HERE == "") HERE <- "."
deg_glob   <- getopt("--deg-glob",  file.path(HERE, "..", "TICTOC_3_factor_model", "Diff_genes_heatmap_*.csv"))
cross_path <- getopt("--crosswalk", file.path(HERE, "..", "crosswalk", "gohir_to_arabidopsis.tsv"))
counts_path<- getopt("--counts",    file.path(HERE, "..", "TICTOC_run1_filteredCounts_v3.csv"))
out_dir    <- getopt("--out",       file.path(HERE, "results"))
padj_cut   <- as.numeric(getopt("--padj", "0.05"))
lfc_cut    <- as.numeric(getopt("--lfc",  "1"))
min_pid    <- as.numeric(getopt("--min-pid", "30"))
use_reg    <- toupper(getopt("--use-regulation", "TRUE")) == "TRUE"  # trust the DEG table's Regulation col
do_kegg    <- toupper(getopt("--kegg", "TRUE")) == "TRUE"            # KEGG needs internet
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(clusterProfiler); library(org.At.tair.db)
has_plots <- requireNamespace("enrichplot", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)

# ---- load crosswalk: Gohir(upper) -> AT locus (one best hit per Gohir) -------
cw <- read.delim(cross_path, stringsAsFactors = FALSE, check.names = TRUE)
if (min_pid > 0 && "pid" %in% names(cw)) cw <- cw[is.na(cw$pid) | cw$pid >= min_pid, ]
g2at <- setNames(cw$arabidopsis_gene, toupper(cw$gohir_gene_upper))
map_to_at <- function(gohir) unique(na.omit(unname(g2at[toupper(gohir)])))
message(sprintf("Crosswalk: %d Gohir->AT pairs (min PID %.0f%%), %d unique AT loci",
                length(g2at), min_pid, length(unique(g2at))))

# ---- universe = AT loci reachable from EXPRESSED cotton genes ----------------
universe <- character(0)
if (file.exists(counts_path)) {
  cm <- read.csv(counts_path, stringsAsFactors = FALSE, check.names = FALSE, nrows = -1)
  universe <- map_to_at(cm[[1]])
  message(sprintf("Universe: %d expressed cotton genes -> %d unique AT loci",
                  nrow(cm), length(universe)))
} else {
  message("WARNING: counts file not found; enriching against whole-genome default universe.")
}

# ---- per-DEG-table enrichment ------------------------------------------------
canon_cols <- function(df) {
  # normalise the (read.csv-mangled) column names to canonical keys
  nm <- tolower(gsub("[^a-z0-9]+", "", tolower(names(df))))
  pick <- function(key) { j <- which(nm == key); if (length(j)) names(df)[j[1]] else NA }
  list(id  = pick("ensemblid"), reg = pick("regulation"),
       lfc = pick("log2foldchange"), padj = pick("adjpval"))
}

run_one <- function(deg_file) {
  tag <- sub("^Diff_genes_heatmap_", "", tools::file_path_sans_ext(basename(deg_file)))
  tag <- gsub("[^A-Za-z0-9._-]+", "_", tag)
  message("\n=== ", tag, " ===")
  df <- read.csv(deg_file, stringsAsFactors = FALSE, check.names = TRUE)
  co <- canon_cols(df)
  if (is.na(co$id)) { message("  no 'Ensembl ID' column; skipping"); return(invisible()) }
  ids <- df[[co$id]]
  padj <- if (!is.na(co$padj)) suppressWarnings(as.numeric(df[[co$padj]])) else rep(NA_real_, nrow(df))
  lfc  <- if (!is.na(co$lfc))  suppressWarnings(as.numeric(df[[co$lfc]]))  else rep(NA_real_, nrow(df))

  # significance + direction
  if (use_reg && !is.na(co$reg)) {
    dir <- toupper(substr(trimws(df[[co$reg]]), 1, 1))          # "U"/"D"
    up   <- ids[dir == "U"]; down <- ids[dir == "D"]
  } else {
    sig  <- !is.na(padj) & padj <= padj_cut & !is.na(lfc) & abs(lfc) >= lfc_cut
    up   <- ids[sig & lfc > 0]; down <- ids[sig & lfc < 0]
  }
  sets <- list(up = up, down = down, all = c(up, down))

  for (dname in names(sets)) {
    at <- map_to_at(sets[[dname]])
    message(sprintf("  %-4s: %d Gohir DEGs -> %d unique AT loci", dname, length(sets[[dname]]), length(at)))
    if (length(at) < 5) { message("    (<5 mapped genes; skipping enrichment)"); next }

    for (ont in c("BP", "MF", "CC")) {
      ego <- tryCatch(enrichGO(gene = at,
                               universe = if (length(universe)) universe else NULL,
                               OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont,
                               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                               readable = TRUE),
                      error = function(e) { message("    GO ", ont, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego),
                  file.path(out_dir, sprintf("GO_%s_%s_%s.csv", ont, dname, tag)), row.names = FALSE)
        if (has_plots) tryCatch({
          p <- enrichplot::dotplot(ego, showCategory = 20) + ggplot2::ggtitle(sprintf("GO:%s %s — %s", ont, dname, tag))
          ggplot2::ggsave(file.path(out_dir, sprintf("GO_%s_%s_%s.pdf", ont, dname, tag)),
                          p, width = 8, height = 8)
        }, error = function(e) NULL)
      }
    }

    if (do_kegg) {
      ek <- tryCatch(enrichKEGG(gene = at, universe = if (length(universe)) universe else NULL,
                                organism = "ath", pvalueCutoff = 0.05),
                     error = function(e) { message("    KEGG failed (needs internet?): ", conditionMessage(e)); NULL })
      if (!is.null(ek) && nrow(as.data.frame(ek)) > 0)
        write.csv(as.data.frame(ek),
                  file.path(out_dir, sprintf("KEGG_%s_%s.csv", dname, tag)), row.names = FALSE)
    }
  }
}

deg_files <- Sys.glob(deg_glob)
if (!length(deg_files)) stop("No DEG files matched: ", deg_glob)
message(sprintf("Found %d DEG table(s).", length(deg_files)))
invisible(lapply(deg_files, run_one))
message("\nDone. Enrichment tables + dotplots in: ", normalizePath(out_dir))
