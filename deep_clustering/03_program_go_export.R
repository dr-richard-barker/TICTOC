#!/usr/bin/env Rscript
# Stage 3 — GO per expression program + export the PhysioSpace reference & AT->Entrez map.
suppressWarnings(suppressMessages({ library(clusterProfiler); library(org.At.tair.db) }))
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 error = function(e) "."); if (is.na(HERE) || HERE == "") HERE <- "."
OUT <- file.path(HERE, "results"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cw <- read.delim(file.path(HERE, "..", "crosswalk", "gohir_to_arabidopsis.tsv"))
g2at <- setNames(cw$arabidopsis_gene, toupper(cw$gohir_gene_upper))
pa <- read.csv(file.path(OUT, "program_assignments.csv"))
bg <- unique(na.omit(g2at[toupper(pa$gene)]))

# GO (BP) per program -> top terms
res <- list()
for (p in sort(unique(pa$program))) {
  genes <- pa$gene[pa$program == p]
  at <- unique(na.omit(g2at[toupper(genes)]))
  ego <- tryCatch(enrichGO(at, universe = bg, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP",
                  pvalueCutoff = 0.1, qvalueCutoff = 0.3, readable = TRUE), error = function(e) NULL)
  df <- if (!is.null(ego)) as.data.frame(ego) else NULL
  top <- if (!is.null(df) && nrow(df)) paste(head(df$Description, 3), collapse = "; ") else "(none)"
  res[[length(res) + 1]] <- data.frame(program = p, n_at = length(at), top_GO = top,
    top1 = if (!is.null(df) && nrow(df)) df$Description[1] else "",
    top1_padj = if (!is.null(df) && nrow(df)) signif(df$p.adjust[1], 3) else NA)
  cat(sprintf("program %2d (%d AT genes): %s\n", p, length(at), top))
}
write.csv(do.call(rbind, res), file.path(OUT, "program_GO.csv"), row.names = FALSE)

# AT locus -> Entrez map for the crosswalk AT genes (for Python stress scoring)
uat <- unique(cw$arabidopsis_gene)
ez <- suppressWarnings(mapIds(org.At.tair.db, keys = uat, column = "ENTREZID", keytype = "TAIR", multiVals = "first"))
write.csv(data.frame(arabidopsis_gene = uat, entrez = unname(ez)), file.path(OUT, "at_to_entrez.csv"), row.names = FALSE)

# Export AT_Stress_Space (Entrez x stress axes) to CSV
if (requireNamespace("PlantPhysioSpace", quietly = TRUE)) {
  e <- new.env(); utils::data("AT_Stress_Space", package = "PlantPhysioSpace", envir = e)
  sp <- as.matrix(get("AT_Stress_Space", envir = e))
  write.csv(data.frame(entrez = rownames(sp), sp, check.names = FALSE),
            file.path(OUT, "AT_Stress_Space.csv"), row.names = FALSE)
  cat("exported AT_Stress_Space:", nrow(sp), "genes x", ncol(sp), "axes\n")
}
cat("stage 3 done\n")
