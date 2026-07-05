#!/usr/bin/env Rscript
# Annotate each WGCNA root module with GO (BP/MF/CC) to derive a human-readable name.
# Reads wgcna/results/module_assignments.csv; writes per-module GO tables + a top-terms summary.
suppressWarnings(suppressMessages({ library(clusterProfiler); library(org.At.tair.db) }))
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 error = function(e) "."); if (is.na(HERE) || HERE == "") HERE <- "."
res <- file.path(HERE, "results")
cw <- read.delim(file.path(HERE, "..", "crosswalk", "gohir_to_arabidopsis.tsv"))
g2at <- setNames(cw$arabidopsis_gene, toupper(cw$gohir_gene_upper))
ma <- read.csv(file.path(res, "module_assignments.csv"))
bg <- unique(na.omit(g2at[toupper(ma$gene)]))

summary_lines <- list()
for (mod in sort(unique(ma$module))) {
  genes <- ma$gene[ma$module == mod]
  at <- unique(na.omit(g2at[toupper(genes)]))
  tops <- character(0)
  for (ont in c("BP", "MF", "CC")) {
    ego <- tryCatch(enrichGO(at, universe = bg, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont,
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE), error = function(e) NULL)
    df <- if (!is.null(ego)) as.data.frame(ego) else NULL
    if (!is.null(df) && nrow(df) > 0) {
      write.csv(df, file.path(res, sprintf("GO_%s_module_%s.csv", ont, mod)), row.names = FALSE)
      tops <- c(tops, paste0(ont, ": ", paste(head(df$Description, 3), collapse = "; ")))
    }
  }
  line <- sprintf("[%s] n=%d, mapped=%d\n    %s", mod, length(genes), length(at),
                  if (length(tops)) paste(tops, collapse = "\n    ") else "(no significant GO enrichment)")
  summary_lines[[mod]] <- line
  cat(line, "\n")
}
writeLines(unlist(summary_lines), file.path(res, "module_GO_top_terms.txt"))
cat("\nwrote module_GO_top_terms.txt\n")
