#!/usr/bin/env Rscript
# =============================================================================
# TICTOC — WGCNA co-expression modules in ROOT, linked to Flight/genotype/RSML traits
# =============================================================================
# Rebuilds the co-expression side of the iDEP v3 report reproducibly, and adds the
# integrative step the manuscript needs (Fig 4 / R4): correlate module eigengenes
# with Treatment, genotype, AND root morphometric traits (from morphometrics/).
#
# Root only (24 samples): tissue dominates PC1 (v4), so co-expression is run within
# the tissue where both the spaceflight response and the RSML traits live.
#
# Outputs (wgcna/results/):
#   module_assignments.csv         gene -> module
#   module_trait_correlation.csv   ME vs Treatment/Genotype/RSML traits (r + p)
#   module_trait_heatmap.pdf
#   flight_module_GO.csv           GO of the top Flight-correlated module
#
# Usage: Rscript run_wgcna_root.R
# Deps: WGCNA, DESeq2, clusterProfiler, org.At.tair.db.
# =============================================================================
suppressWarnings(suppressMessages({ library(WGCNA); library(DESeq2) }))
options(stringsAsFactors = FALSE); enableWGCNAThreads <- try(WGCNA::enableWGCNAThreads(), silent = TRUE)
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 error = function(e) "."); if (is.na(HERE) || HERE == "") HERE <- "."
out_dir <- file.path(HERE, "results"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
counts_path <- file.path(HERE, "..", "TICTOC_run1_filteredCounts_v3.csv")
cross_path  <- file.path(HERE, "..", "crosswalk", "gohir_to_arabidopsis.tsv")
summ_path   <- file.path(HERE, "..", "morphometrics", "rsml_traits_summary.csv")

# ---- root VST, top-variable genes -------------------------------------------
cm <- read.csv(counts_path, row.names = 1, check.names = FALSE)
cts <- as.matrix(cm); mode(cts) <- "integer"
grp <- sub("(__r[0-9]+|[0-9]+)$", "", colnames(cts)); colnames(cts) <- make.unique(colnames(cts), sep = "__r")
geno <- sub("(Root|Shoot).*", "", grp); tis <- sub(".*?(Root|Shoot).*", "\\1", grp)
trt <- sub(".*(Flight|Ground)$", "\\1", grp)
root <- tis == "Root"
cd <- data.frame(row.names = colnames(cts)[root], genotype = geno[root], treatment = trt[root])
dds <- DESeqDataSetFromMatrix(cts[, root], cd, ~1); v <- assay(vst(dds, blind = TRUE))
topg <- head(order(-rowVars(v)), 5000)
datExpr <- t(v[topg, ])
cat(sprintf("root samples=%d, genes=%d\n", nrow(datExpr), ncol(datExpr)))

# ---- soft power + modules ---------------------------------------------------
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
sp <- sft$powerEstimate; if (is.na(sp)) sp <- 12
cat("soft power =", sp, "\n")
cor <- WGCNA::cor      # avoid the WGCNA/other-package cor() masking that breaks blockwiseModules
net <- blockwiseModules(datExpr, power = sp, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
                        maxBlockSize = 6000, verbose = 0)
cor <- stats::cor      # restore base cor for the trait-correlation step below
moduleColors <- labels2colors(net$colors)
write.csv(data.frame(gene = colnames(datExpr), module = moduleColors),
          file.path(out_dir, "module_assignments.csv"), row.names = FALSE)
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
cat("modules:", paste(table(moduleColors), names(table(moduleColors)), collapse = " | "), "\n")

# ---- trait matrix: Treatment, genotype, RSML root traits (group-level) ------
traits <- data.frame(Flight = as.integer(cd$treatment == "Flight"),
                     AVPOX  = as.integer(cd$genotype %in% c("A68", "D130")),
                     A68    = as.integer(cd$genotype == "A68"),
                     D130   = as.integer(cd$genotype == "D130"),
                     row.names = rownames(cd))
if (file.exists(summ_path)) {                              # group-level RSML traits (day 6)
  s <- read.csv(summ_path); s6 <- s[s$day == 6, ]
  key <- paste(s6$genotype, ifelse(s6$condition == "FL", "Flight", "Ground"))
  for (tr in c("total_length_native_mean", "primary_length_native_mean", "n_lateral_mean")) {
    lut <- setNames(s6[[tr]], key)
    traits[[sub("_mean$", "", tr)]] <- lut[paste(cd$genotype, cd$treatment)]
  }
}

mt_cor <- cor(MEs, traits, use = "pairwise.complete.obs")
mt_p   <- corPvalueStudent(mt_cor, nrow(datExpr))
res <- data.frame(module = rownames(mt_cor))
for (cl in colnames(mt_cor)) { res[[paste0("r_", cl)]] <- round(mt_cor[, cl], 3); res[[paste0("p_", cl)]] <- signif(mt_p[, cl], 3) }
write.csv(res, file.path(out_dir, "module_trait_correlation.csv"), row.names = FALSE)

pdf(file.path(out_dir, "module_trait_heatmap.pdf"), width = 8, height = 7)
labeledHeatmap(Matrix = mt_cor, xLabels = colnames(traits), yLabels = rownames(mt_cor),
               ySymbols = rownames(mt_cor), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = paste0(round(mt_cor, 2), "\n(", signif(mt_p, 1), ")"),
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1, 1), main = "Root module–trait correlations")
dev.off()

# ---- GO of the top Flight-correlated module ---------------------------------
topmod <- rownames(mt_cor)[which.max(abs(mt_cor[, "Flight"]))]
mod_col <- sub("^ME", "", topmod)
mod_genes <- colnames(datExpr)[moduleColors == mod_col]
cat("top Flight module:", topmod, "r=", round(mt_cor[topmod, "Flight"], 2), "(", length(mod_genes), "genes)\n")
ok <- requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.At.tair.db", quietly = TRUE)
if (ok && file.exists(cross_path)) {
  cw <- read.delim(cross_path); g2at <- setNames(cw$arabidopsis_gene, toupper(cw$gohir_gene_upper))
  at <- unique(na.omit(g2at[toupper(mod_genes)]))
  bg <- unique(na.omit(g2at[toupper(colnames(datExpr))]))
  ego <- try(clusterProfiler::enrichGO(at, universe = bg, OrgDb = org.At.tair.db::org.At.tair.db,
             keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE), silent = TRUE)
  if (!inherits(ego, "try-error") && !is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file.path(out_dir, sprintf("flight_module_%s_GO.csv", mod_col)), row.names = FALSE)
    cat("  GO top:", paste(head(as.data.frame(ego)$Description, 4), collapse = " | "), "\n")
  }
}
cat("WGCNA done ->", out_dir, "\n")
