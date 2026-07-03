#!/usr/bin/env Rscript
# =============================================================================
# TICTOC — PhysioSpace stress-pattern decoding of cotton spaceflight contrasts
# =============================================================================
# Mirrors the OSD-767 tomato pipeline (roadmap 4.5): project cotton Flight-vs-Ground
# fold changes onto Arabidopsis stress reference spaces to get PhysioScores per
# genotype x tissue axis.
#
# Method (as used for OSD-767; Lenz et al. 2013; Hadizadeh Esfahani et al. 2021):
#   1. Build an INPUT matrix = per-(Genotype x Tissue) Flight - Ground fold change,
#      computed from VST-transformed counts (genes x contrasts, in Gohir space).
#   2. Re-index to Arabidopsis: Gohir -> AT locus (crosswalk) -> AT Entrez
#      (org.At.tair.db), collapsing many-Gohir -> one-Entrez by mean.
#   3. calculatePhysioMap(Input, Space, GenesRatio=0.05, TTEST=FALSE,
#      ImputationMethod="PCA") against each reference space.
#   4. Write PhysioScores + a heatmap per reference space.
#
# YOU MUST SUPPLY THE REFERENCE SPACES (not in this repo — same files as OSD-767):
#   AT_Stress_Space, AT_Stress_Space_Meta, AT_Stress_Space_RNASeq
#   as .rds matrices (rownames = Arabidopsis Entrez IDs, cols = stress axes).
#   Point --ref-dir at the folder holding them (see physiospace/README.md).
#
# Usage:
#   Rscript run_physiospace.R --ref-dir /path/to/reference_spaces \
#       [--counts ../TICTOC_run1_filteredCounts_v3.csv] \
#       [--crosswalk ../crosswalk/gohir_to_arabidopsis.tsv] \
#       [--out results] [--min-pid 30] [--genes-ratio 0.05]
#
# Install (once):
#   install.packages(c("BiocManager","remotes","pheatmap"))
#   BiocManager::install(c("DESeq2","org.At.tair.db"))
#   remotes::install_github("JRC-COMBINE/PhysioSpaceMethods")
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[[i + 1]]
}
HERE <- tryCatch(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), error = function(e) ".")
if (is.na(HERE) || HERE == "") HERE <- "."
counts_path <- getopt("--counts",    file.path(HERE, "..", "TICTOC_run1_filteredCounts_v3.csv"))
cross_path  <- getopt("--crosswalk", file.path(HERE, "..", "crosswalk", "gohir_to_arabidopsis.tsv"))
ref_dir     <- getopt("--ref-dir",   NA_character_)
out_dir     <- getopt("--out",       file.path(HERE, "results"))
min_pid     <- as.numeric(getopt("--min-pid", "30"))
genes_ratio <- as.numeric(getopt("--genes-ratio", "0.05"))
if (is.na(ref_dir)) stop("Provide --ref-dir pointing at the AT stress reference-space .rds files.")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

need <- c("DESeq2", "org.At.tair.db", "PhysioSpaceMethods")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "),
                       "\nSee the install block in this script's header.")
suppressmessages <- function(x) suppressWarnings(suppressMessages(x))
suppressmessages(library(DESeq2)); suppressmessages(library(org.At.tair.db)); suppressmessages(library(PhysioSpaceMethods))

# ---- 1. counts -> sample table (parse group labels in the header) -----------
cm <- read.csv(counts_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
cts <- as.matrix(cm); mode(cts) <- "integer"
grp <- sub("[0-9]+$", "", colnames(cts))                       # A68RootFlight1 -> A68RootFlight
parse1 <- function(g) {
  geno <- sub("(Root|Shoot).*", "", g)
  tis  <- sub(".*?(Root|Shoot).*", "\\1", g)
  trt  <- sub(".*(Flight|Ground)$", "\\1", g)
  c(genotype = geno, tissue = tis, treatment = trt)
}
meta <- as.data.frame(t(vapply(grp, parse1, character(3))), stringsAsFactors = FALSE)
rownames(meta) <- colnames(cts)
meta[] <- lapply(meta, factor)
message(sprintf("Counts: %d genes x %d samples; groups: %s",
                nrow(cts), ncol(cts), paste(levels(interaction(meta$genotype, meta$tissue)), collapse = ", ")))

# ---- 2. VST, then Flight - Ground mean per (genotype x tissue) ---------------
dds <- DESeqDataSetFromMatrix(cts, meta, design = ~ 1)
vsd <- assay(vst(dds, blind = TRUE))                           # genes x samples (Gohir rows)
combos <- unique(meta[, c("genotype", "tissue")])
fc <- sapply(seq_len(nrow(combos)), function(i) {
  g <- combos$genotype[i]; t <- combos$tissue[i]
  fl <- rownames(meta)[meta$genotype == g & meta$tissue == t & meta$treatment == "Flight"]
  gr <- rownames(meta)[meta$genotype == g & meta$tissue == t & meta$treatment == "Ground"]
  if (!length(fl) || !length(gr)) return(rep(NA_real_, nrow(vsd)))
  rowMeans(vsd[, fl, drop = FALSE]) - rowMeans(vsd[, gr, drop = FALSE])
})
colnames(fc) <- paste0(combos$genotype, "_", combos$tissue, "_FlightVsGround")
fc <- fc[, colSums(!is.na(fc)) > 0, drop = FALSE]
message("Input contrasts: ", paste(colnames(fc), collapse = ", "))

# ---- 3. Gohir -> AT locus -> AT Entrez, collapse many->one by mean ----------
cw <- read.delim(cross_path, stringsAsFactors = FALSE)
if (min_pid > 0 && "pid" %in% names(cw)) cw <- cw[is.na(cw$pid) | cw$pid >= min_pid, ]
g2at <- setNames(cw$arabidopsis_gene, toupper(cw$gohir_gene_upper))
at_of <- g2at[toupper(rownames(fc))]                           # Gohir rowname -> AT locus (may be NA)
uat <- unique(unname(at_of[!is.na(at_of)]))
at2ez <- suppressWarnings(mapIds(org.At.tair.db, keys = uat,
                                 column = "ENTREZID", keytype = "TAIR", multiVals = "first"))
map_entrez <- setNames(at2ez[at_of], rownames(fc))             # Gohir rowname -> Entrez (may be NA)
keep <- !is.na(map_entrez)
fc2 <- fc[keep, , drop = FALSE]; ez <- as.character(map_entrez[keep])
sums <- rowsum(fc2, group = ez, na.rm = TRUE)                  # sum rows sharing an Entrez
cnts <- as.vector(table(ez)[rownames(sums)])                  # group sizes, aligned to sums' rownames
input <- sums / cnts                                           # MEAN across homoeologs/paralogs
message(sprintf("Re-indexed to Arabidopsis Entrez: %d genes x %d contrasts",
                nrow(input), ncol(input)))

# ---- 4. PhysioScores against each reference space ---------------------------
ref_files <- list.files(ref_dir, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
if (!length(ref_files)) stop("No .rds reference spaces found in ", ref_dir)
for (rf in ref_files) {
  space <- readRDS(rf); sp_name <- tools::file_path_sans_ext(basename(rf))
  message("\n== Reference space: ", sp_name, " (", nrow(space), " genes x ", ncol(space), " axes) ==")
  shared <- intersect(rownames(input), rownames(space))
  message("  shared Entrez genes with input: ", length(shared))
  scores <- tryCatch(
    calculatePhysioMap(InputData = input, Space = space,
                       GenesRatio = genes_ratio, TTEST = FALSE, ImputationMethod = "PCA"),
    error = function(e) { message("  calculatePhysioMap failed: ", conditionMessage(e)); NULL })
  if (is.null(scores)) next
  write.csv(scores, file.path(out_dir, sprintf("PhysioScores_%s.csv", sp_name)))
  if (requireNamespace("pheatmap", quietly = TRUE)) tryCatch(
    pheatmap::pheatmap(scores, main = paste("PhysioScores —", sp_name), cluster_cols = FALSE,
                       filename = file.path(out_dir, sprintf("PhysioScores_%s.pdf", sp_name))),
    error = function(e) NULL)
}
message("\nDone. PhysioScores + heatmaps in: ", normalizePath(out_dir))
