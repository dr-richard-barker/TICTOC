#!/usr/bin/env Rscript
# Named module-trait heatmap (Fig 4) + module_names.csv, from GO-derived human-readable names.
suppressWarnings(suppressMessages({ library(WGCNA); library(DESeq2) }))
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 error = function(e) "."); if (is.na(HERE) || HERE == "") HERE <- "."

names_lut <- c(turquoise = "Signalling & isoprenoid metab.", blue = "Translation & ribosome",
               brown = "Defence & ubiquitin signalling", green = "Photosynthesis (light rxn)",
               yellow = "Metal transport & phenylpropanoid", grey = "Unassigned")
write.csv(data.frame(module = names(names_lut), name = names_lut),
          file.path(HERE, "..", "wgcna", "results", "module_names.csv"), row.names = FALSE)

cm <- read.csv(file.path(HERE, "..", "TICTOC_run1_filteredCounts_v3.csv"), row.names = 1, check.names = FALSE)
cts <- as.matrix(cm); mode(cts) <- "integer"
grp <- sub("(__r[0-9]+|[0-9]+)$", "", colnames(cts)); colnames(cts) <- make.unique(colnames(cts), sep = "__r")
geno <- sub("(Root|Shoot).*", "", grp); tis <- sub(".*?(Root|Shoot).*", "\\1", grp)
trt <- sub(".*(Flight|Ground)$", "\\1", grp); root <- tis == "Root"
cd <- data.frame(row.names = colnames(cts)[root], genotype = geno[root], treatment = trt[root])
v <- assay(vst(DESeqDataSetFromMatrix(cts[, root], cd, ~1), blind = TRUE))
ma <- read.csv(file.path(HERE, "..", "wgcna", "results", "module_assignments.csv"))
v <- v[rownames(v) %in% ma$gene, ]; datExpr <- t(v)
modcol <- setNames(ma$module, ma$gene)[colnames(datExpr)]
MEs <- moduleEigengenes(datExpr, modcol)$eigengenes
grpkey <- paste(cd$genotype, cd$treatment); ME_g <- apply(MEs, 2, function(c) tapply(c, grpkey, mean))
s <- read.csv(file.path(HERE, "..", "morphometrics", "rsml_traits_summary.csv")); s6 <- s[s$day == 6, ]
tkey <- paste(s6$genotype, ifelse(s6$condition == "FL", "Flight", "Ground"))
tr <- data.frame(row.names = rownames(ME_g),
  Flight = as.integer(grepl("Flight", rownames(ME_g))), AVPOX = as.integer(grepl("A68|D130", rownames(ME_g))),
  `Total root length` = setNames(s6$total_length_native_mean, tkey)[rownames(ME_g)],
  `Primary length` = setNames(s6$primary_length_native_mean, tkey)[rownames(ME_g)],
  `Lateral count` = setNames(s6$n_lateral_mean, tkey)[rownames(ME_g)], check.names = FALSE)

mt <- cor(ME_g, tr, use = "pairwise.complete.obs"); mtp <- corPvalueStudent(mt, nrow(ME_g))
lab <- paste0(sub("^ME", "", rownames(mt)), " — ", names_lut[sub("^ME", "", rownames(mt))])
ord <- order(-mt[, "Flight"])
draw_hm <- function() {
  par(mar = c(7, 13, 3, 1))
  labeledHeatmap(Matrix = mt[ord, ], xLabels = colnames(tr), yLabels = lab[ord],
                 colors = blueWhiteRed(50), textMatrix = paste0(round(mt[ord, ], 2)),
                 setStdMargins = FALSE, cex.text = 0.7, cex.lab = 0.8, zlim = c(-1, 1),
                 main = "Root co-expression modules vs Flight, AVP-OX and root architecture (group-level, n=6)")
}
pdf(file.path(HERE, "..", "manuscript", "Fig4_module_trait_named.pdf"), width = 8.5, height = 5.5); draw_hm(); dev.off()
png(file.path(HERE, "..", "manuscript", "Fig4_module_trait_named.png"), width = 8.5, height = 5.5, units = "in", res = 140); draw_hm(); dev.off()
cat("wrote Fig4_module_trait_named.pdf + module_names.csv\n")
