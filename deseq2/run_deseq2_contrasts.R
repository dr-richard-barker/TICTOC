#!/usr/bin/env Rscript
# =============================================================================
# TICTOC — full DESeq2 contrast set from the filtered count matrix
# =============================================================================
# Fits one model over all 48 samples (12 groups = Genotype x Tissue x Treatment,
# 4 reps) and extracts the canonical contrasts, writing DEG tables in the same
# column format as TICTOC_3_factor_model/Diff_genes_heatmap_*.csv so they flow
# straight into go_analysis/. Also writes a DEG-count summary (Table 2 / R3).
#
# Contrasts:
#   Spaceflight (Flight - Ground) within each Genotype x Tissue        (6)
#   AVP-OX interaction: (line FL-GC) - (WT FL-GC), per tissue          (4)   [A68,D130 x Root,Shoot]
#
# Usage:
#   Rscript run_deseq2_contrasts.R [--counts ../TICTOC_run1_filteredCounts_v3.csv]
#                                  [--out contrasts] [--padj 0.05] [--lfc 1]
# Deps: DESeq2 (+ ashr if available for LFC shrinkage; falls back to MLE LFC).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[[i + 1]] }
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 error = function(e) ".")
if (is.na(HERE) || HERE == "") HERE <- "."
counts_path <- getopt("--counts", file.path(HERE, "..", "TICTOC_run1_filteredCounts_v3.csv"))
out_dir     <- getopt("--out", file.path(HERE, "contrasts"))
padj_cut    <- as.numeric(getopt("--padj", "0.05"))
lfc_cut     <- as.numeric(getopt("--lfc", "1"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
suppressWarnings(suppressMessages(library(DESeq2)))
has_ashr <- requireNamespace("ashr", quietly = TRUE)

# ---- counts -> coldata (parse group labels; replicate cols share a label) ---
cm  <- read.csv(counts_path, row.names = 1, check.names = FALSE)
cts <- as.matrix(cm); mode(cts) <- "integer"
grp <- sub("(__r[0-9]+|[0-9]+)$", "", colnames(cts))
colnames(cts) <- make.unique(colnames(cts), sep = "__r")
geno <- sub("(Root|Shoot).*", "", grp)
tis  <- sub(".*?(Root|Shoot).*", "\\1", grp)
trt  <- sub(".*(Flight|Ground)$", "\\1", grp)
group <- factor(paste(geno, tis, trt, sep = "_"))
coldata <- data.frame(row.names = colnames(cts), group = group)
cat(sprintf("counts %d x %d; %d groups: %s\n", nrow(cts), ncol(cts), nlevels(group),
            paste(levels(group), collapse = ", ")))

dds <- DESeqDataSetFromMatrix(cts, coldata, design = ~ 0 + group)
dds <- DESeq(dds)
rn <- resultsNames(dds)                                    # e.g. "groupWT_Root_Ground"
gvec <- function(lbl) { v <- setNames(numeric(length(rn)), rn); v[paste0("group", lbl)] <- NA; v }
mk <- function(pos, neg) {                                 # numeric contrast: +pos groups, -neg groups
  v <- setNames(numeric(length(rn)), rn)
  for (p in pos) v[paste0("group", p)] <- v[paste0("group", p)] + 1
  for (n in neg) v[paste0("group", n)] <- v[paste0("group", n)] - 1
  v
}

get_res <- function(con) {
  if (has_ashr) tryCatch(lfcShrink(dds, contrast = con, type = "ashr", quiet = TRUE),
                         error = function(e) results(dds, contrast = con))
  else results(dds, contrast = con)
}

write_deg <- function(res, tag) {
  df <- as.data.frame(res); df <- df[!is.na(df$padj), ]
  reg <- ifelse(df$padj < padj_cut & df$log2FoldChange >=  lfc_cut, "Up",
         ifelse(df$padj < padj_cut & df$log2FoldChange <= -lfc_cut, "Down", "NS"))
  out <- data.frame(Regulation = reg, `Ensembl ID` = rownames(df),
                    `log2 Fold Change` = round(df$log2FoldChange, 4), `Adj.Pval` = df$padj,
                    baseMean = round(df$baseMean, 2), pvalue = df$pvalue, check.names = FALSE)
  out <- out[order(out$Adj.Pval), ]
  write.csv(out, file.path(out_dir, sprintf("Diff_genes_heatmap_%s.csv", tag)), row.names = FALSE)
  c(up = sum(reg == "Up"), down = sum(reg == "Down"))
}

genos <- c("WT", "A68", "D130"); tissues <- c("Root", "Shoot")
summary_rows <- list()
# 6 simple spaceflight contrasts
for (g in genos) for (t in tissues) {
  tag <- sprintf("%s_%s_FlightVsGround", g, t)
  cnt <- write_deg(get_res(mk(paste(g, t, "Flight", sep = "_"), paste(g, t, "Ground", sep = "_"))), tag)
  summary_rows[[tag]] <- c(contrast = tag, type = "spaceflight", cnt)
  cat(sprintf("  %-28s up=%d down=%d\n", tag, cnt["up"], cnt["down"]))
}
# 4 AVP-OX interaction contrasts
for (g in c("A68", "D130")) for (t in tissues) {
  tag <- sprintf("%svsWT_%s_FlightInteraction", g, t)
  con <- mk(pos = c(paste(g, t, "Flight", sep = "_"), paste("WT", t, "Ground", sep = "_")),
            neg = c(paste(g, t, "Ground", sep = "_"), paste("WT", t, "Flight", sep = "_")))
  cnt <- write_deg(get_res(con), tag)
  summary_rows[[tag]] <- c(contrast = tag, type = "interaction", cnt)
  cat(sprintf("  %-28s up=%d down=%d\n", tag, cnt["up"], cnt["down"]))
}

summ <- do.call(rbind, lapply(summary_rows, function(r) as.data.frame(as.list(r), stringsAsFactors = FALSE)))
write.csv(summ, file.path(out_dir, "DEG_counts_summary.csv"), row.names = FALSE)
cat("shrinkage:", if (has_ashr) "ashr" else "none (MLE LFC)",
    "| wrote", length(summary_rows), "contrasts +>", file.path(out_dir, "DEG_counts_summary.csv"), "\n")
