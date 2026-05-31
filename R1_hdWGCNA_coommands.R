#!/usr/bin/env Rscript
# hdWGCNA + TF-target analysis workflow
# Usage: Rscript run_dataset.R <dataset_number>
# e.g.: Rscript run_dataset.R 6

# ============================================================
# LIBRARY PATH - add system site library (contains Genomic* packages)
# ============================================================
.libPaths(c(.libPaths(), "/usr/local/lib/R/site-library"))

# Progress logger: writes to stderr (always unbuffered, unlike stdout when piped)
log_progress <- function(...) cat(..., file = stderr())

# ============================================================
# SEURAT 5 COMPATIBILITY PATCH FOR hdWGCNA
# hdWGCNA 0.4.8 calls GetAssayData(slot=) which is defunct in Seurat 5.
# Override the S3 method dispatch via registerS3method() to silently
# convert 'slot' to 'layer', making the call work transparently.
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})
if(packageVersion("SeuratObject") >= "5.0.0") {
  message("Applying Seurat 5 compatibility patch for hdWGCNA...")
  .compat_GetAssayData_Seurat <- function(object, assay = NULL, slot = NULL,
                                          layer = NULL, ...) {
    # Silently convert deprecated 'slot' to 'layer'
    if (!is.null(slot) && !inherits(slot, "deprecated") && is.null(layer)) {
      layer <- slot
    }
    if (is.null(layer)) layer <- "data"
    object <- SeuratObject::UpdateSlots(object = object)
    assay  <- assay %||% SeuratObject::DefaultAssay(object = object)
    assay_obj <- object[[assay]]
    SeuratObject::LayerData(assay_obj, layer = layer)
  }
  registerS3method("GetAssayData", "Seurat", .compat_GetAssayData_Seurat,
                   envir = getNamespace("SeuratObject"))
  message("Seurat 5 compatibility patch applied.")
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) stop("Usage: Rscript run_dataset.R <dataset_number>")
dataset_num <- as.integer(args[1])

# ============================================================
# DATASET-SPECIFIC CONFIGURATION
# ============================================================
BASE_DIR <- "/tank/projects/maxim_kovalev/TFtarget"
# Set working directory to BASE_DIR so hdWGCNA tom_outdir is relative and paths don't double
setwd(BASE_DIR)

config <- list(
  # Macaca fascicularis: use human genome (hg38) + EnsDb as proxy.
  # JASPAR vertebrate motifs and gene names are highly conserved between
  # macaque and human microglia. Species-specific BSgenome (NCBI.6.0) tried
  # if installed, falls back to hg38 if not.
  `6` = list(
    h5ad       = file.path(BASE_DIR, "Dataset_6.h5ad"),
    outdir     = file.path(BASE_DIR, "Dataset_6"),
    species    = "Macaca fascicularis",
    genome_pkg = "BSgenome.Mfascicularis.NCBI.6.0",
    genome_pkg_fallback = "BSgenome.Hsapiens.UCSC.hg38",
    genome_str = "Mfas6.0",
    genome_str_fallback = "hg38",
    ensdb_query = c("EnsDb", "Macaca fascicularis"),
    ensdb_pkg_fallback = "EnsDb.Hsapiens.v86",
    jaspar_tax  = "vertebrates",
    sample_col  = "sample"
  ),
  # Macaca mulatta: same proxy strategy
  `7` = list(
    h5ad       = file.path(BASE_DIR, "Dataset_7.h5ad"),
    outdir     = file.path(BASE_DIR, "Dataset_7"),
    species    = "Macaca mulatta",
    genome_pkg = "BSgenome.Mmulatta.UCSC.rheMac10",
    genome_pkg_fallback = "BSgenome.Hsapiens.UCSC.hg38",
    genome_str = "rheMac10",
    genome_str_fallback = "hg38",
    ensdb_query = c("EnsDb", "Macaca mulatta"),
    ensdb_pkg_fallback = "EnsDb.Hsapiens.v86",
    jaspar_tax  = "vertebrates",
    sample_col  = "sample"
  ),
  `8` = list(
    h5ad    = file.path(BASE_DIR, "Dataset_8.h5ad"),
    outdir  = file.path(BASE_DIR, "Dataset_8"),
    species = "Homo sapiens",
    genome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
    genome_str = "hg38",
    ensdb_pkg  = "EnsDb.Hsapiens.v86",
    jaspar_tax  = "vertebrates",
    sample_col  = "sample"
  ),
  `9` = list(
    h5ad    = file.path(BASE_DIR, "Dataset_9.h5ad"),
    outdir  = file.path(BASE_DIR, "Dataset_9"),
    species = "Homo sapiens",
    genome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
    genome_str = "hg38",
    ensdb_pkg  = "EnsDb.Hsapiens.v86",
    jaspar_tax  = "vertebrates",
    sample_col  = "sample_batch"
  )
)

cfg <- config[[as.character(dataset_num)]]
if(is.null(cfg)) stop(paste("Unknown dataset:", dataset_num))

outdir   <- cfg$outdir
plotdir  <- file.path(outdir, "plots")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
# Switch working directory to the dataset's own folder so WGCNA consensusTOM
# temporary block files (consensusTOM-block.N.RData) are written there and don't
# conflict between parallel dataset runs.
setwd(outdir)

log_progress("\n", strrep("=", 60), "\n")
log_progress("Dataset:", dataset_num, "-", cfg$species, "\n")
log_progress("Output:", outdir, "\n")
log_progress(strrep("=", 60), "\n\n")

# ============================================================
# LOAD LIBRARIES
# ============================================================
log_progress("[1/7] Loading libraries...\n")
suppressPackageStartupMessages({
  library(hdWGCNA)
  library(Seurat)
  library(Matrix)
  library(anndata)
  library(JASPAR2024)
  library(TFBSTools)
  library(motifmatchr)
  library(GenomicRanges)
  library(ensembldb)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Load BSgenome (with fallback to human genome for macaque datasets)
genome_loaded <- FALSE
if(!is.null(cfg$genome_pkg) && cfg$genome_pkg %in% rownames(installed.packages())) {
  tryCatch({
    library(cfg$genome_pkg, character.only = TRUE)
    genome_obj <- get(cfg$genome_pkg)
    genome_str_use <- cfg$genome_str
    log_progress("  Genome:", cfg$genome_pkg, "\n")
    genome_loaded <- TRUE
  }, error = function(e) log_progress("  WARNING: Failed to load", cfg$genome_pkg, ":", e$message, "\n"))
}
if(!genome_loaded) {
  fallback_pkg <- if(!is.null(cfg$genome_pkg_fallback)) cfg$genome_pkg_fallback else "BSgenome.Hsapiens.UCSC.hg38"
  fallback_str <- if(!is.null(cfg$genome_str_fallback)) cfg$genome_str_fallback else "hg38"
  log_progress("  Using fallback genome:", fallback_pkg, "\n")
  library(fallback_pkg, character.only = TRUE)
  genome_obj <- get(fallback_pkg)
  genome_str_use <- fallback_str
}

# Load or retrieve EnsDb (with fallback to human EnsDb for macaque datasets)
EnsDb_loaded <- FALSE
if(!is.null(cfg$ensdb_pkg)) {
  tryCatch({
    library(cfg$ensdb_pkg, character.only = TRUE)
    EnsDb_obj <- get(cfg$ensdb_pkg)
    log_progress("  EnsDb:", cfg$ensdb_pkg, "\n")
    EnsDb_loaded <- TRUE
  }, error = function(e) log_progress("  WARNING: Failed to load", cfg$ensdb_pkg, ":", e$message, "\n"))
}
if(!EnsDb_loaded && !is.null(cfg$ensdb_query)) {
  log_progress("  Trying AnnotationHub for", cfg$species, "...\n")
  tryCatch({
    library(AnnotationHub)
    ah <- AnnotationHub(localHub = FALSE)
    ah_res <- query(ah, cfg$ensdb_query)
    if(length(ah_res) > 0) {
      latest_id <- tail(names(ah_res), 1)
      EnsDb_obj  <- ah[[latest_id]]
      log_progress("  Using EnsDb:", latest_id, "\n")
      EnsDb_loaded <- TRUE
    }
  }, error = function(e) log_progress("  AnnotationHub failed:", e$message, "\n"))
}
if(!EnsDb_loaded) {
  fallback_ensdb <- if(!is.null(cfg$ensdb_pkg_fallback)) cfg$ensdb_pkg_fallback else "EnsDb.Hsapiens.v86"
  log_progress("  Using fallback EnsDb:", fallback_ensdb, "\n")
  library(fallback_ensdb, character.only = TRUE)
  EnsDb_obj <- get(fallback_ensdb)
  # Update genome_str_use to match the fallback EnsDb
  if(fallback_ensdb == "EnsDb.Hsapiens.v86") genome_str_use <- "hg38"
}

log_progress("Libraries loaded.\n\n")

# ============================================================
# STEP 1: CONVERT H5AD TO SEURAT OBJECT (with checkpoint)
# ============================================================
raw_rds <- file.path(outdir, "seurat_raw.rds")
if(file.exists(raw_rds)) {
  log_progress("[2/7] Loading existing seurat_raw.rds (checkpoint)...\n")
  seurat_obj <- readRDS(raw_rds)
  log_progress("  Loaded:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")
} else {
  log_progress("[2/7] Converting h5ad to Seurat...\n")

  adata <- read_h5ad(cfg$h5ad)
  log_progress("  Loaded:", nrow(adata$obs), "cells x", nrow(adata$var), "genes\n")

  # Get matrices — anndata returns (cells x genes), Seurat needs (genes x cells)
  # X = log1p normalized, counts layer = raw counts
  log1p_mat  <- t(adata$X)                        # genes x cells
  counts_mat <- t(adata$layers[["counts"]])        # genes x cells

  # Make gene names unique (Seurat requires unique rownames)
  if(any(duplicated(rownames(counts_mat)))){
    log_progress("  WARNING: Duplicate gene names found, making unique\n")
    rownames(log1p_mat)  <- make.unique(rownames(log1p_mat))
    rownames(counts_mat) <- make.unique(rownames(counts_mat))
  }

  # Make cell names unique.
  # IMPORTANT: when obs index has duplicates, py_to_r sets rownames(meta) as "1","2","3"...
  # (not actual barcodes). Check colnames(counts_mat) for duplicates instead.
  meta <- as.data.frame(adata$obs)
  cell_names_raw <- colnames(counts_mat)
  if(any(duplicated(cell_names_raw))){
    log_progress("  WARNING: Duplicate cell names found (", sum(duplicated(cell_names_raw)), "), making unique\n")
    unique_cells <- make.unique(cell_names_raw)
    rownames(meta)       <- unique_cells
    colnames(log1p_mat)  <- unique_cells
    colnames(counts_mat) <- unique_cells
  } else {
    rownames(meta) <- cell_names_raw
  }
  cell_names <- colnames(counts_mat)

  # Create Seurat object (RNA assay with raw counts)
  seurat_obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta)

  # Set log1p normalized data in the data layer (Seurat 5 API)
  # Seurat may have modified gene names (e.g., underscores -> dashes).
  # Row order is preserved from counts_mat, so align by position not by name.
  rownames(log1p_mat) <- rownames(seurat_obj)
  LayerData(seurat_obj, "data") <- log1p_mat[, colnames(seurat_obj)]

  # Add pre-computed dimensionality reductions from h5ad
  pca_emb  <- adata$obsm[["X_pca"]];  rownames(pca_emb)  <- cell_names
  umap_emb <- adata$obsm[["X_umap"]]; rownames(umap_emb) <- cell_names
  colnames(pca_emb)  <- paste0("PC_",   seq_len(ncol(pca_emb)))
  colnames(umap_emb) <- paste0("UMAP_", seq_len(ncol(umap_emb)))
  seurat_obj[["pca"]]  <- CreateDimReducObject(pca_emb,  key = "PC_",   assay = "RNA")
  seurat_obj[["umap"]] <- CreateDimReducObject(umap_emb, key = "UMAP_", assay = "RNA")
  # NOTE: ScaleData is NOT called — PCA/UMAP are preloaded, and WGCNA uses log-normalized data.

  log_progress("  Seurat object:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")
  log_progress("  Layers:", paste(SeuratObject::Layers(seurat_obj[["RNA"]]), collapse=", "), "\n")
  log_progress("  Samples:", paste(unique(seurat_obj@meta.data[[cfg$sample_col]]), collapse=", "), "\n\n")

  # Save checkpoint
  saveRDS(seurat_obj, raw_rds)
  log_progress("  seurat_raw.rds saved.\n")
}

wgcna_rds <- file.path(outdir, "seurat_wgcna.rds")
if(file.exists(wgcna_rds)) {
  log_progress("[3-5/7] Loading existing seurat_wgcna.rds (WGCNA checkpoint)...\n")
  seurat_obj <- readRDS(wgcna_rds)
  log_progress("  WGCNA checkpoint loaded.\n\n")
} else {

# ============================================================
# STEP 2: hdWGCNA SETUP & METACELLS
# ============================================================
log_progress("[3/7] hdWGCNA: SetupForWGCNA, metacells...\n")

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction    = 0.001,
  wgcna_name  = "wgcna"
)
log_progress("  Genes selected for WGCNA:", length(GetWGCNAGenes(seurat_obj)), "\n")

# Metacells — k=25 is standard; for small datasets use smaller k
n_cells <- ncol(seurat_obj)
k_val <- ifelse(n_cells < 2000, 10, 25)
min_cells_val <- max(100, round(n_cells / 50))

log_progress("  Constructing metacells: k =", k_val, ", min_cells =", min_cells_val, "\n")
seurat_obj <- MetacellsByGroups(
  seurat_obj  = seurat_obj,
  group.by    = c(cfg$sample_col),
  k           = k_val,
  reduction   = "pca",
  ident.group = cfg$sample_col,
  min_cells   = min_cells_val
)
seurat_obj <- NormalizeMetacells(seurat_obj)

# ============================================================
# STEP 3: SET DATA EXPRESSION & TEST SOFT POWERS
# ============================================================
log_progress("[4/7] hdWGCNA: SetDatExpr, TestSoftPowers...\n")

# Use only groups that survived MetacellsByGroups filtering
metacell_obj <- GetMetacellObject(seurat_obj)
valid_samples <- unique(metacell_obj@meta.data[[cfg$sample_col]])
log_progress("  Samples with sufficient cells for metacells:", paste(valid_samples, collapse=", "), "\n")

seurat_obj <- SetDatExpr(
  seurat_obj,
  group.by   = cfg$sample_col,
  group_name = valid_samples,
  assay      = "RNA",
  slot       = "data"
)

seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")

power_table <- GetPowerTable(seurat_obj)
write.csv(power_table, file.path(outdir, "soft_power_table.csv"), row.names = FALSE)

pdf(file.path(plotdir, "soft_power.pdf"), width = 12, height = 8)
print(PlotSoftPowers(seurat_obj))
dev.off()

# Select soft power: first power where SFT.R.sq >= 0.8
sp_candidates <- power_table$Power[power_table$SFT.R.sq >= 0.80]
soft_power <- if(length(sp_candidates) > 0) min(sp_candidates) else {
  diffs <- diff(power_table$SFT.R.sq)
  power_table$Power[which.min(diffs) + 1]
}
soft_power <- max(1, min(soft_power, 30))
log_progress("  Selected soft power:", soft_power, "\n")
log_progress("  (SFT.R.sq at selected power:", power_table$SFT.R.sq[power_table$Power == soft_power], ")\n")

# ============================================================
# STEP 4: CONSTRUCT NETWORK & MODULES
# ============================================================
log_progress("[5/7] hdWGCNA: ConstructNetwork, Modules...\n")

seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power  = soft_power,
  setDatExpr  = FALSE,
  tom_outdir  = ".",  # relative to getwd()=outdir; each dataset has its own working dir
  overwrite_tom = TRUE
)

tryCatch({
  pdf(file.path(plotdir, "dendrogram.pdf"), width = 10, height = 6)
  PlotDendrogram(seurat_obj, main = paste0("Dataset_", dataset_num, " hdWGCNA Dendrogram"))
  dev.off()
}, error = function(e) { dev.off(); log_progress("  Dendrogram plot skipped:", e$message, "\n") })

seurat_obj <- ModuleEigengenes(seurat_obj)  # no Harmony (ScaleData not required; TF network uses metacell expr)
seurat_obj <- ModuleConnectivity(seurat_obj)
seurat_obj <- RunModuleUMAP(seurat_obj, n_hubs = 10, n_neighbors = 15, min_dist = 0.1)

modules <- GetModules(seurat_obj)
n_modules <- length(unique(modules$module[modules$module != "grey"]))
log_progress("  Number of non-grey modules:", n_modules, "\n")
write.csv(modules, file.path(outdir, "modules.csv"), row.names = FALSE)

tryCatch({
  pdf(file.path(plotdir, "module_umap.pdf"), width = 10, height = 8)
  ModuleUMAPPlot(seurat_obj, label_hubs = 5)
  dev.off()
}, error = function(e) {
  dev.off()
  log_progress("  ModuleUMAPPlot skipped:", e$message, "\n")
})

tryCatch({
  MEs <- GetMEs(seurat_obj, harmonized = FALSE)
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  module_cols <- colnames(MEs)
  if(length(module_cols) > 0){
    pdf(file.path(plotdir, "module_feature_plots.pdf"), width = 12, height = 8)
    for(i in seq(1, min(length(module_cols), 20), by = 6)){
      cols_chunk <- module_cols[i:min(i+5, length(module_cols))]
      p <- FeaturePlot(seurat_obj, features = cols_chunk, reduction = "umap", ncol = 3, order = TRUE) &
           theme(plot.title = element_text(size=10))
      print(p)
    }
    dev.off()
  }
}, error = function(e) {
  log_progress("  Feature plots skipped:", e$message, "\n")
})

hub_genes <- GetHubGenes(seurat_obj, n_hubs = 25)
write.csv(hub_genes, file.path(outdir, "hub_genes.csv"), row.names = FALSE)
log_progress("  Hub genes saved.\n")

saveRDS(seurat_obj, wgcna_rds)
log_progress("  seurat_wgcna.rds saved.\n\n")

} # end WGCNA block

# ============================================================
# STEP 5: TF MOTIF SCANNING
# ============================================================
log_progress("[6/7] TF analysis: MotifScan...\n")

# Get JASPAR2024 motifs (vertebrate CORE collection)
# JASPAR2024 package changed API: JASPAR2024() returns S4 object with @db SQLite path
# TFBSTools 1.44.0 getMatrixSet() works with SQLiteConnection but not JASPAR2024 class
# So we use the SQLite path directly.
jaspar_db_obj <- JASPAR2024::JASPAR2024()
jaspar_sqlite_path <- jaspar_db_obj@db
jaspar_conn <- DBI::dbConnect(RSQLite::SQLite(), jaspar_sqlite_path)
pfm_set <- getMatrixSet(
  jaspar_conn,
  opts = list(
    collection  = "CORE",
    tax_group   = cfg$jaspar_tax,
    all_versions = FALSE
  )
)
DBI::dbDisconnect(jaspar_conn)
log_progress("  JASPAR2024 motifs loaded:", length(pfm_set), "motifs\n")

log_progress("  Using genome:", genome_str_use, "with EnsDb:", class(EnsDb_obj), "\n")
# Run MotifScan
seurat_obj <- MotifScan(
  seurat_obj,
  pfm            = pfm_set,
  EnsDb          = EnsDb_obj,
  species_genome = genome_str_use,
  wgcna_name     = "wgcna"
)

# Add gene_name column to motif table (required by ConstructTFNetwork)
motif_df <- GetMotifs(seurat_obj)
# JASPAR motif names are typically TF gene names (e.g., "SP1", "CTCF")
# Clean up: some motifs have "(var.X)" suffix or "::" for composite motifs
motif_df$gene_name <- gsub("\\(var\\..*\\)$", "", motif_df$motif_name)
motif_df$gene_name <- gsub("::.*$", "", motif_df$gene_name)  # take first TF of composite
motif_df$gene_name <- trimws(motif_df$gene_name)
# Use uppercase for consistency
motif_df$gene_name <- toupper(motif_df$gene_name)
seurat_obj <- SetMotifs(seurat_obj, motif_df)

# Overlap motifs with modules
seurat_obj <- OverlapModulesMotifs(seurat_obj)

# Save motif overlap results
motif_overlap <- GetMotifOverlap(seurat_obj)
write.csv(motif_overlap, file.path(outdir, "motif_module_overlap.csv"), row.names = FALSE)

# Plot motif enrichment
pdf(file.path(plotdir, "motif_enrichment.pdf"), width = 12, height = 8)
tryCatch(
  MotifEnrichmentPlot(seurat_obj, seurat_ref = NULL, plot_size = c(5,12)),
  error = function(e) {
    log_progress("  MotifEnrichmentPlot failed:", e$message, "\n")
    grid::grid.text("MotifEnrichmentPlot not available")
  }
)
dev.off()

log_progress("  MotifScan complete.\n\n")

# ============================================================
# STEP 6: CONSTRUCT TF NETWORK
# ============================================================
log_progress("[7/7] TF analysis: ConstructTFNetwork...\n")

source("/tank/projects/maxim_kovalev/TFtarget/ConstructTFNetwork.R")

model_params <- list(
  booster    = "gbtree",
  objective  = "reg:squarederror",
  eval_metric = "rmse",
  eta        = 0.1,
  max_depth  = 6,
  subsample  = 0.7,
  colsample_bytree = 0.7,
  nthread    = 8
)

seurat_obj <- ConstructTFNetwork(
  seurat_obj,
  model_params = model_params,
  nfold        = 5,
  wgcna_name   = "wgcna"
)

# Save TF network results
tf_network  <- GetTFNetwork(seurat_obj)
tf_eval     <- GetTFEval(seurat_obj)
write.csv(tf_network, file.path(outdir, "tf_network.csv"),  row.names = FALSE)
write.csv(tf_eval,    file.path(outdir, "tf_eval.csv"),     row.names = FALSE)

# Top TF-target pairs per module
if(!is.null(tf_network) && nrow(tf_network) > 0){
  top_tfs <- tf_network %>%
    group_by(gene) %>%
    slice_max(Gain, n = 5) %>%
    ungroup()
  write.csv(top_tfs, file.path(outdir, "top_tf_targets.csv"), row.names = FALSE)
  log_progress("  Top TF-target pairs saved.\n")
}

# Plot TF network (hub TFs)
pdf(file.path(plotdir, "tf_network.pdf"), width = 12, height = 10)
tryCatch(
  PlotTFNetwork(seurat_obj, n_hubs = 10, vertex.label.cex = 0.7),
  error = function(e) {
    log_progress("  PlotTFNetwork failed:", e$message, "\n")
    grid::grid.text("TF Network plot not available")
  }
)
dev.off()

# ============================================================
# FINAL SAVE
# ============================================================
log_progress("\nSaving final Seurat object...\n")
saveRDS(seurat_obj, file.path(outdir, "seurat_final.rds"))

log_progress("\n", strrep("=", 60), "\n")
log_progress("Dataset_", dataset_num, " COMPLETE!\n")
log_progress("Output files in:", outdir, "\n")
log_progress(strrep("=", 60), "\n\n")

# Print summary
log_progress("Summary:\n")
log_progress("  Cells:", ncol(seurat_obj), "\n")
log_progress("  Genes:", nrow(seurat_obj), "\n")
log_progress("  WGCNA genes:", length(GetWGCNAGenes(seurat_obj)), "\n")
mods <- GetModules(seurat_obj)
n_mods <- length(unique(mods$module[mods$module != "grey"]))
log_progress("  Modules (non-grey):", n_mods, "\n")
log_progress("  Soft power used:", soft_power, "\n")
log_progress("  Genome:", genome_str_use, "\n")
log_progress("  Species:", cfg$species, "\n")
if(!is.null(tf_network)) log_progress("  TF-target pairs:", nrow(tf_network), "\n")
