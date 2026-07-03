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
# REFERENCE SPACES come from the PlantPhysioSpace data package (Entrez-indexed AT stress spaces).
#   Default = the three used by OSD-767: AT_Stress_Space, AT_Stress_Space_Meta, AT_Stress_Space_RNASeq.
#   Override with --spaces "A,B,C" (any data() object in PlantPhysioSpace), or point --ref-dir at a
#   folder of .rds matrices if you have local copies.
#
# Usage:
#   Rscript run_physiospace.R
#       [--spaces AT_Stress_Space,AT_Stress_Space_Meta,AT_Stress_Space_RNASeq]
#       [--ref-dir /path/to/local/rds]        # optional: bypass the data package
#       [--counts ../TICTOC_run1_filteredCounts_v3.csv]
#       [--crosswalk ../crosswalk/gohir_to_arabidopsis.tsv]
#       [--out results] [--min-pid 30] [--genes-ratio 0.05]
#
# Refs: Lenz et al. 2013 PLoS ONE (10.1371/journal.pone.0077627);
#       Hadizadeh Esfahani et al. 2021 Plant Physiol. 187:1795 (Plant PhysioSpace).
#
# Install (once):
#   install.packages(c("BiocManager","remotes"))
#   BiocManager::install(c("DESeq2","org.At.tair.db"))
#   remotes::install_github("JRC-COMBINE/PhysioSpaceMethods")   # method
#   remotes::install_github("JRC-COMBINE/PlantPhysioSpace")     # data (AT/rice/soy/wheat spaces)
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
ref_dir     <- getopt("--ref-dir",   NA_character_)   # optional: local .rds override
spaces_arg  <- getopt("--spaces",    "AT_Stress_Space,AT_Stress_Space_Meta,AT_Stress_Space_RNASeq")
out_dir     <- getopt("--out",       file.path(HERE, "results"))
min_pid     <- as.numeric(getopt("--min-pid", "30"))
genes_ratio <- as.numeric(getopt("--genes-ratio", "0.05"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

need <- c("DESeq2", "org.At.tair.db", "PhysioSpaceMethods")
if (is.na(ref_dir)) need <- c(need, "PlantPhysioSpace")   # only needed when loading spaces from the package
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

# ---- 4. resolve reference spaces (PlantPhysioSpace package, or local .rds) ---
space_names <- trimws(strsplit(spaces_arg, ",")[[1]])
load_space <- function(nm) {
  if (!is.na(ref_dir)) {                                   # local override
    f <- file.path(ref_dir, paste0(nm, ".rds"))
    if (!file.exists(f)) { message("  (no local file ", f, ")"); return(NULL) }
    return(readRDS(f))
  }
  e <- new.env()
  ok <- tryCatch({ utils::data(list = nm, package = "PlantPhysioSpace", envir = e); TRUE },
                 error = function(err) FALSE)
  if (!ok || !exists(nm, envir = e)) { message("  (space '", nm, "' not found in PlantPhysioSpace)"); return(NULL) }
  get(nm, envir = e)
}

# ---- 5. PhysioScores against each reference space ---------------------------
for (sp_name in space_names) {
  space <- load_space(sp_name); if (is.null(space)) next
  space <- as.matrix(space)
  message("\n== Reference space: ", sp_name, " (", nrow(space), " genes x ", ncol(space), " axes) ==")
  message("  shared Entrez genes with input: ", length(intersect(rownames(input), rownames(space))))
  scores <- tryCatch(
    calculatePhysioMap(InputData = input, Space = space,
                       GenesRatio = genes_ratio, TTEST = FALSE, ImputationMethod = "PCA"),
    error = function(e) { message("  calculatePhysioMap failed: ", conditionMessage(e)); NULL })
  if (is.null(scores)) next
  write.csv(scores, file.path(out_dir, sprintf("PhysioScores_%s.csv", sp_name)))
  # prefer the package's own heatmap; fall back to pheatmap
  plotted <- FALSE
  if ("PhysioHeatmap" %in% getNamespaceExports("PhysioSpaceMethods")) tryCatch({
    pdf(file.path(out_dir, sprintf("PhysioScores_%s.pdf", sp_name)), width = 9, height = 7)
    PhysioSpaceMethods::PhysioHeatmap(PhysioResults = scores, main = paste("PhysioScores —", sp_name),
                                      SymmetricColoring = TRUE, SpaceClustering = FALSE)
    dev.off(); plotted <- TRUE
  }, error = function(e) { try(dev.off(), silent = TRUE); message("  PhysioHeatmap failed: ", conditionMessage(e)) })
  if (!plotted && requireNamespace("pheatmap", quietly = TRUE)) tryCatch(
    pheatmap::pheatmap(scores, main = paste("PhysioScores —", sp_name), cluster_cols = FALSE,
                       filename = file.path(out_dir, sprintf("PhysioScores_%s.pdf", sp_name))),
    error = function(e) NULL)
}
message("\nDone. PhysioScores + heatmaps in: ", normalizePath(out_dir))
