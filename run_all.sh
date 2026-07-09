#!/usr/bin/env bash
# TICTOC — reproduce the whole analysis end to end. See REPRODUCE.md for prerequisites.
# Usage:  RSCRIPT=/path/to/Rscript  bash run_all.sh    (RSCRIPT defaults to `Rscript` on PATH)
set -euo pipefail
cd "$(dirname "$0")"
RS="${RSCRIPT:-Rscript}"
PY="${PYTHON:-python}"
log(){ printf '\n=== %s ===\n' "$*"; }

log "1/9 ortholog crosswalk";        "$PY" crosswalk/build_gohir_to_arabidopsis.py
log "2/9 morphometrics";             "$PY" morphometrics/extract_rsml_traits.py
                                     "$PY" morphometrics/summarise_traits.py
                                     "$PY" morphometrics/morphometric_stats.py
log "3/9 DESeq2 contrasts";          "$RS" deseq2/run_deseq2_contrasts.R
log "4/9 DESeq2 factorial (V4)";     "$RS" deseq2/run_factorial_v4.R
log "5/9 GO / KEGG";                 "$RS" go_analysis/run_go_clusterprofiler.R \
                                        --deg-glob "deseq2/contrasts/Diff_genes_heatmap_*.csv" \
                                        --out go_analysis/results_full
log "6/9 PhysioSpace";               "$RS" physiospace/run_physiospace.R --static TRUE \
                                        --out physiospace/results_static
log "7/9 WGCNA";                     "$RS" wgcna/run_wgcna_root.R
                                     "$RS" wgcna/annotate_modules_go.R
log "8/9 integration";              "$RS" integration/integrate_root_expression.R
                                     "$RS" integration/make_named_heatmap.R
log "9/10 figures";                  "$RS" deseq2/make_pca.R
                                     "$PY" manuscript/figures/make_figures.py
log "10/10 deep-learning programs (R8, Fig 8-9)"
                                     "$PY" deep_clustering/01_autoencoder_programs.py
                                     "$PY" deep_clustering/02_celltype_markers.py
                                     "$RS" deep_clustering/03_program_go_export.R
                                     "$PY" deep_clustering/04_program_stress.py
                                     "$PY" deep_clustering/05_figures.py
log "DONE — outputs under each stage's results/ folder; figures in manuscript/figures/"
