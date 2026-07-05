#!/usr/bin/env Rscript
# =============================================================================
# TICTOC V4 — reproducible factorial interaction model (improves on iDEP v3)
# =============================================================================
# The original TICTOC_markdown_v3.html (iDEP 0.95) fit a Treatment x Genotype x
# Tissue model with all two-way interactions (refs Ground/WT/Root). This script
# reproduces that model in plain DESeq2 (versioned, no hard-coded paths) and
# quantifies each interaction term, plus PC-factor correlations — the two most
# useful v3 results to carry into the manuscript.
#
# Outputs (deseq2/v4/):
#   interaction_DEG_counts.csv  — significant DEGs per interaction term
#   pc_factor_correlation.csv   — which factor each PC tracks (Kruskal-Wallis)
#   v4_summary.txt
#
# Usage: Rscript run_factorial_v4.R [--counts ../TICTOC_run1_filteredCounts_v3.csv]
# Deps: DESeq2.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[[i + 1]] }
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 error = function(e) ".")
if (is.na(HERE) || HERE == "") HERE <- "."
counts_path <- getopt("--counts", file.path(HERE, "..", "TICTOC_run1_filteredCounts_v3.csv"))
out_dir <- getopt("--out", file.path(HERE, "v4")); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
padj_cut <- as.numeric(getopt("--padj", "0.05")); lfc_cut <- as.numeric(getopt("--lfc", "1"))
suppressWarnings(suppressMessages(library(DESeq2)))

cm <- read.csv(counts_path, row.names = 1, check.names = FALSE)
cts <- as.matrix(cm); mode(cts) <- "integer"
grp <- sub("(__r[0-9]+|[0-9]+)$", "", colnames(cts)); colnames(cts) <- make.unique(colnames(cts), sep = "__r")
cd <- data.frame(row.names = colnames(cts),
                 Treatment = relevel(factor(sub(".*(Flight|Ground)$", "\\1", grp)), "Ground"),
                 Genotype  = relevel(factor(sub("(Root|Shoot).*", "", grp)), "WT"),
                 Tissue    = relevel(factor(sub(".*?(Root|Shoot).*", "\\1", grp)), "Root"))

# depth-vs-factor QC (the v3 warning)
depth <- colSums(cts)
qc <- sapply(c("Treatment","Genotype","Tissue"), function(f) kruskal.test(depth, cd[[f]])$p.value)

dds <- DESeqDataSetFromMatrix(cts, cd, design = ~ Treatment * Genotype * Tissue)
dds <- DESeq(dds)
rn <- resultsNames(dds)
int_terms <- grep("\\.", rn, value = TRUE)              # interaction coefficients contain a "."
cnt <- t(sapply(int_terms, function(t) {
  r <- results(dds, name = t)
  c(up = sum(!is.na(r$padj) & r$padj < padj_cut & r$log2FoldChange >= lfc_cut),
    down = sum(!is.na(r$padj) & r$padj < padj_cut & r$log2FoldChange <= -lfc_cut))
}))
ic <- data.frame(term = rownames(cnt), cnt, total = rowSums(cnt), row.names = NULL)
ic <- ic[order(-ic$total), ]
write.csv(ic, file.path(out_dir, "interaction_DEG_counts.csv"), row.names = FALSE)

# PC-factor correlation (which factor each PC tracks)
vsd <- vst(dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)[head(order(-rowVars(assay(vsd))), 2000), ]))
pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pcf <- do.call(rbind, lapply(1:5, function(i) {
  data.frame(PC = paste0("PC", i), pct_var = pv[i],
             p_Treatment = kruskal.test(pca$x[, i], cd$Treatment)$p.value,
             p_Genotype  = kruskal.test(pca$x[, i], cd$Genotype)$p.value,
             p_Tissue    = kruskal.test(pca$x[, i], cd$Tissue)$p.value)
}))
write.csv(pcf, file.path(out_dir, "pc_factor_correlation.csv"), row.names = FALSE)

sink(file.path(out_dir, "v4_summary.txt"))
cat("TICTOC V4 factorial model — ~ Treatment*Genotype*Tissue (refs Ground/WT/Root)\n\n")
cat("Depth-vs-factor QC (Kruskal p; low p = read-depth confound):\n"); print(round(qc, 4))
cat("\nInteraction-term DEG counts (padj<0.05, |LFC|>=1):\n"); print(ic, row.names = FALSE)
cat("\nPC-factor correlation (Kruskal p):\n"); print(pcf, row.names = FALSE)
sink()
cat("V4 done ->", out_dir, "\n"); print(ic, row.names = FALSE); print(pcf[, c("PC","pct_var","p_Treatment","p_Tissue")], row.names = FALSE)
