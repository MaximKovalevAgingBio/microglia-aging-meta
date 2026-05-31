#!/usr/bin/env Rscript
# ==============================================================================
# 02: Genome-Wide Single-Cell Co-Expression Estimation via CS-CORE
# This script processes a filtered Seurat object to compute an all-against-all
# co-expression matrix adjusted for sequencing depth and technical noise.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CSCORE)
})

# Pipeline Configuration (Standardized for Dataset 1)
ds  <- 1
rds <- "/tank/projects/maxim_kovalev/TFtarget/Dataset_1/seurat_final.rds"
out <- "/tank/projects/maxim_kovalev/CS-CORE"

cat(sprintf("[DS%d] Loading seurat object...\n", ds))
obj <- readRDS(rds)
cat(sprintf("[DS%d] %d genes x %d cells  (assay: %s)\n",
            ds, nrow(obj), ncol(obj), DefaultAssay(obj)))

all_genes <- rownames(obj)
cat(sprintf("[DS%d] Extracting count matrix for %d genes...\n", ds, length(all_genes)))

# Extract counts and sequence depth for C++ optimization engine
count_matrix <- t(as.matrix(LayerData(obj, assay = "RNA", layer = "counts")))
seq_depth    <- obj$nCount_RNA

cat(sprintf("[DS%d] Running CSCORE_IRLS_cpp on %d genes x %d cells...\n",
            ds, ncol(count_matrix), nrow(count_matrix)))

result <- CSCORE_IRLS_cpp(
  count_matrix, seq_depth,
  IRLS_par = list(n_iter = 10, eps = 0.05, verbose = TRUE, conv = "max")
)

# Name rows and columns with gene names
rownames(result$est)     <- all_genes
colnames(result$est)     <- all_genes
rownames(result$p_value) <- all_genes
colnames(result$p_value) <- all_genes

# ==============================================================================
# OPTIMIZED FDR MULTIPLE TEST CORRECTION (Benjamini-Hochberg)
# ==============================================================================
cat(sprintf("[DS%d] Applying FDR multiple test correction across the matrix...\n", ds))

# Extract upper triangle to save RAM (since the matrix is symmetric)
upper_tri_indices <- upper.tri(result$p_value, diag = FALSE)
raw_p_vector      <- result$p_value[upper_tri_indices]

# Compute Benjamini-Hochberg adjusted p-values
adj_p_vector      <- p.adjust(raw_p_vector, method = "BH")

# Create a clean matrix for adjusted p-values and populate it symmetrically
fdr_matrix <- matrix(NA, nrow = nrow(result$p_value), ncol = ncol(result$p_value))
fdr_matrix[upper_tri_indices] <- adj_p_vector

# Reflect upper triangle to lower triangle for perfect square symmetry
fdr_matrix <- t(fdr_matrix)
fdr_matrix[upper_tri_indices] <- adj_p_vector

# Set diagonal elements to 0 (self-coexpression is always highly significant)
diag(fdr_matrix) <- 0

# Restore metadata names
rownames(fdr_matrix) <- all_genes
colnames(fdr_matrix) <- all_genes

# ==============================================================================
# WRITING MATRICES TO DISK
# ==============================================================================
cat(sprintf("[DS%d] Writing output CSV matrices to disk...\n", ds))

# Save estimates (co-expression coefficients)
write.csv(result$est,
          file.path(out, sprintf("Dataset%d_CSCORE_FULL_matrix_est.csv", ds)),
          row.names = TRUE)

# Save the newly calculated FDR-adjusted p-values matrix
write.csv(fdr_matrix,
          file.path(out, sprintf("Dataset%d_CSCORE_FULL_matrix_p.csv", ds)),
          row.names = TRUE)

cat(sprintf("[DS%d] Done. Pipeline complete. Output shape: %d x %d\n",
            ds, nrow(result$est), ncol(result$est)))