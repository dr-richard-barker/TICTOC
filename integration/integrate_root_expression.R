#!/usr/bin/env Rscript
# =============================================================================
# TICTOC — integrate root transcriptome (WGCNA modules) with root architecture
# =============================================================================
# Goal: identify genes / modules / pathways whose expression co-varies with the
# measured root-growth changes, as a correlational (hypothesis-generating) model.
#
# STATISTICAL FRAMING (read first):
#   RNA-seq (4 replicate libraries per Genotype x Tissue x Treatment group) is NOT
#   individually paired to the imaged plants (RSML is per-plant). The only valid
#   join is at the GROUP level: 6 root genotype x treatment groups. All expression-
#   trait correlations below are therefore computed on n = 6 group means (not 24
#   pseudo-replicated samples) and are HYPOTHESIS-GENERATING, not causal. A causal
#   integrated model (sPLS/DIABLO, or mediation Treatment->expression->trait) would
#   require individually-paired expression+imaging or more groups/timepoints.
#
# Layers:
#   1. Module eigengene <-> trait correlation at group level (n=6).
#   2. Module Membership (MM) within-modality (24 samples) -> hub genes.
#   3. Gene Significance (GS) = group-level cor(gene, root trait) -> growth-correlated genes.
#   4. Hub-and-growth genes = high |MM| in a growth module AND high |GS|.
#
# Outputs (integration/results/): module_trait_grouplevel.csv, gene_GS_MM.csv,
#   hub_growth_genes.csv, integration_summary.csv
# Deps: WGCNA, DESeq2.
# =============================================================================
suppressWarnings(suppressMessages({ library(WGCNA); library(DESeq2) }))
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 error = function(e) "."); if (is.na(HERE) || HERE == "") HERE <- "."
out <- file.path(HERE, "results"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
GROWTH <- "total_length_native"                             # headline root-growth trait

cm <- read.csv(file.path(HERE, "..", "TICTOC_run1_filteredCounts_v3.csv"), row.names = 1, check.names = FALSE)
cts <- as.matrix(cm); mode(cts) <- "integer"
grp <- sub("(__r[0-9]+|[0-9]+)$", "", colnames(cts)); colnames(cts) <- make.unique(colnames(cts), sep = "__r")
geno <- sub("(Root|Shoot).*", "", grp); tis <- sub(".*?(Root|Shoot).*", "\\1", grp)
trt <- sub(".*(Flight|Ground)$", "\\1", grp); root <- tis == "Root"
cd <- data.frame(row.names = colnames(cts)[root], genotype = geno[root], treatment = trt[root])
v <- assay(vst(DESeqDataSetFromMatrix(cts[, root], cd, ~1), blind = TRUE))

ma <- read.csv(file.path(HERE, "..", "wgcna", "results", "module_assignments.csv"))
v <- v[rownames(v) %in% ma$gene, ]; datExpr <- t(v)         # samples x module genes
modcol <- setNames(ma$module, ma$gene)[colnames(datExpr)]
MEs <- orderMEs(moduleEigengenes(datExpr, modcol)$eigengenes)

# group means (n = 6 genotype x treatment groups)
grpkey <- paste(cd$genotype, cd$treatment)
gmean <- function(M) apply(M, 2, function(col) tapply(col, grpkey, mean))
ME_g <- gmean(MEs); expr_g <- gmean(datExpr)

# group-level traits (day-6 RSML means, root)
s <- read.csv(file.path(HERE, "..", "morphometrics", "rsml_traits_summary.csv")); s6 <- s[s$day == 6, ]
tkey <- paste(s6$genotype, ifelse(s6$condition == "FL", "Flight", "Ground"))
traits_g <- data.frame(row.names = rownames(ME_g),
  Flight = as.integer(grepl("Flight", rownames(ME_g))),
  AVPOX  = as.integer(grepl("A68|D130", rownames(ME_g))))
for (tr in c("total_length_native", "primary_length_native", "n_lateral"))
  traits_g[[tr]] <- setNames(s6[[paste0(tr, "_mean")]], tkey)[rownames(ME_g)]

# 1. module-trait at group level (n=6)
mt <- cor(ME_g, traits_g, use = "pairwise.complete.obs")
mtp <- corPvalueStudent(mt, nrow(ME_g))
mtab <- data.frame(module = sub("^ME", "", rownames(mt)))
for (cl in colnames(mt)) { mtab[[paste0("r_", cl)]] <- round(mt[, cl], 3); mtab[[paste0("p_", cl)]] <- signif(mtp[, cl], 2) }
write.csv(mtab, file.path(out, "module_trait_grouplevel.csv"), row.names = FALSE)

# 2+3. MM (24 samples, within-modality) and GS (group-level, n=6)
MM <- cor(datExpr, MEs, use = "pairwise.complete.obs")
GS <- cor(expr_g, traits_g[[GROWTH]], use = "pairwise.complete.obs")
gm <- data.frame(gene = colnames(datExpr), module = modcol,
                 MM_own = mapply(function(g, m) MM[g, paste0("ME", m)], colnames(datExpr), modcol),
                 GS_growth = GS[, 1])
gm <- gm[order(-abs(gm$GS_growth)), ]
write.csv(gm, file.path(out, "gene_GS_MM.csv"), row.names = FALSE)

# 4. hub-and-growth genes in growth-associated modules
growth_mods <- sub("^ME", "", rownames(mt))[abs(mt[, GROWTH]) > 0.5]
hub <- subset(gm, module %in% growth_mods & abs(MM_own) > 0.7 & abs(GS_growth) > 0.7)
hub <- hub[order(hub$module, -abs(hub$GS_growth)), ]
write.csv(hub, file.path(out, "hub_growth_genes.csv"), row.names = FALSE)

# 5. integration summary (per module)
isum <- mtab[, c("module", "r_Flight", "r_AVPOX", "r_total_length_native", "r_n_lateral")]
isum$n_genes <- as.integer(table(modcol)[isum$module])
isum$n_hub_growth <- as.integer(sapply(isum$module, function(m) sum(hub$module == m)))
isum <- isum[order(-abs(isum$r_total_length_native)), ]
write.csv(isum, file.path(out, "integration_summary.csv"), row.names = FALSE)

cat("growth-associated modules (|r_len|>0.5):", paste(growth_mods, collapse = ", "), "\n")
cat("hub-and-growth genes:", nrow(hub), "in", paste(unique(hub$module), collapse = ", "), "\n")
print(isum, row.names = FALSE)
cat("\ntop growth-correlated genes:\n"); print(head(gm[, c("gene","module","MM_own","GS_growth")], 10), row.names = FALSE)
