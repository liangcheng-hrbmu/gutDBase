# File: scRNAseq_pipeline.R
# Description: Standardized scRNA-seq processing and analysis pipeline for non-batch-corrected data
# Tools / package versions: Seurat v3.2.3, MAST v1.16.0, clusterProfiler v3.18.1
# Note: This pipeline processes individual datasets without batch correction

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)        # v3.2.3 - Single-cell analysis
  library(dplyr)         # Data manipulation
  library(Matrix)        # Sparse matrix operations
  library(ggplot2)       # Data visualization
  library(cowplot)       # Plot arrangement
  library(MAST)          # v1.16.0 - Differential expression testing
  library(clusterProfiler) # Functional enrichment analysis
  library(org.Hs.eg.db)  # Human gene annotation (use org.Mm.eg.db for mouse)
})

# -----------------------------
# USER CONFIGURATION SECTION
# -----------------------------
# Set input and output directories
cellranger_dirs <- list.dirs(path = "./cellranger_outputs", recursive = FALSE)
output_dir <- "./results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Quality control thresholds
min_features_per_cell <- 200    # Minimum features per cell
min_cells_per_gene <- 3         # Minimum cells expressing a gene
max_mito_fraction <- 0.20       # Maximum mitochondrial percentage (20%)

# Differential expression analysis parameters
deg_pval_cutoff <- 0.05         # P-value cutoff for DEGs
deg_log2fc_top_pct <- 0.05      # Top 5% by absolute log2FC for enrichment

# Clustering parameters
n_pcs_to_use <- NULL            # If NULL, use ElbowPlot to determine
resolution <- 0.8               # Clustering resolution

# Species specification (human or mouse)
species <- "human"              # Change to "mouse" if needed

# -----------------------------
# HELPER FUNCTIONS
# -----------------------------
# Calculate mitochondrial percentage
get_mito_pct <- function(seurat_obj, species = c("human","mouse")){
  species <- match.arg(species)
  if(species == "human"){
    mito.pattern <- "^MT-"
  } else {
    mito.pattern <- "^mt-"
  }
  mito.features <- grep(pattern = mito.pattern, 
                       x = rownames(seurat_obj@assays$RNA@counts), 
                       value = TRUE)
  if(length(mito.features) == 0) return(rep(0, ncol(seurat_obj)))
  mito_counts <- Matrix::colSums(seurat_obj@assays$RNA@counts[mito.features, , drop=FALSE])
  pct <- mito_counts / Matrix::colSums(seurat_obj@assays$RNA@counts)
  return(pct)
}

# -----------------------------
# MAIN PROCESSING PIPELINE
# -----------------------------
process_single_dataset <- function(dir, species = "human") {
  sample_name <- basename(dir)
  message("Processing sample: ", sample_name)
  
  # 1. DATA LOADING
  # Locate and load Cell Ranger output
  mat_dir <- file.path(dir, "filtered_feature_bc_matrix")
  if(!dir.exists(mat_dir)){
    mat_dir <- file.path(dir, "outs", "filtered_feature_bc_matrix")
  }
  if(!dir.exists(mat_dir)){
    warning("Cannot find filtered_feature_bc_matrix in ", dir, "; skipping")
    return(NULL)
  }
  
  counts <- Read10X(data.dir = mat_dir)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, 
                                  min.cells = 0, min.features = 0)
  
  # 2. QUALITY CONTROL
  # Filter cells and genes
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= min_features_per_cell)
  seurat_obj <- seurat_obj[Matrix::rowSums(seurat_obj@assays$RNA@counts > 0) >= min_cells_per_gene, ]
  
  # Calculate and filter by mitochondrial percentage
  mito_pct <- get_mito_pct(seurat_obj, species = species)
  seurat_obj$percent.mt <- mito_pct * 100
  seurat_obj <- subset(seurat_obj, subset = percent.mt <= 100 * max_mito_fraction)
  
  # Save QC summary
  qc_summary <- data.frame(
    sample = sample_name,
    cells_after_qc = ncol(seurat_obj),
    genes_after_qc = nrow(seurat_obj),
    mean_features = mean(seurat_obj$nFeature_RNA),
    mean_umis = mean(seurat_obj$nCount_RNA),
    mean_mito_pct = mean(seurat_obj$percent.mt)
  )
  write.csv(qc_summary, file = file.path(output_dir, paste0(sample_name, "_qc_summary.csv")), 
            row.names = FALSE)
  
  # 3. NORMALIZATION
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                             scale.factor = 10000)
  
  # 4. FEATURE SELECTION
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                    nfeatures = 2000)
  
  # 5. SCALING
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  
  # 6. DIMENSIONALITY REDUCTION (PCA)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), 
                      npcs = 50, verbose = FALSE)
  
  # Determine optimal PCs using elbow plot
  if(is.null(n_pcs_to_use)) {
    elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
    ggsave(file.path(output_dir, paste0(sample_name, "_elbow_plot.png")), elbow_plot)
    n_pcs_to_use <- 30  # Default value, adjust based on elbow plot
  }
  
  # 7. CLUSTERING
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:n_pcs_to_use)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  # 8. NON-LINEAR DIMENSIONALITY REDUCTION (UMAP/t-SNE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:n_pcs_to_use)
  seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:n_pcs_to_use)
  
  # 9. CELL TYPE ANNOTATION (placeholder - requires reference data)
  # seurat_obj <- annotate_cell_types(seurat_obj, reference_dataset)
  
  # 10. DIFFERENTIAL EXPRESSION ANALYSIS
  # Find markers for all clusters
  cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, 
                                   logfc.threshold = 0.25, min.pct = 0.25)
  write.csv(cluster_markers, file.path(output_dir, paste0(sample_name, "_cluster_markers.csv")),
            row.names = FALSE)
  
  # 11. VISUALIZATION
  # UMAP colored by cluster
  umap_cluster <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
  ggsave(file.path(output_dir, paste0(sample_name, "_umap_clusters.png")), umap_cluster)
  
  # Feature plots for top markers
  top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
  feature_plot <- FeaturePlot(seurat_obj, features = top_markers$gene[1:4], 
                             reduction = "umap")
  ggsave(file.path(output_dir, paste0(sample_name, "_feature_plot.png")), feature_plot)
  
  # 12. SAVE PROCESSED OBJECT
  saveRDS(seurat_obj, file = file.path(output_dir, paste0(sample_name, "_processed.rds")))
  
  message("Completed processing: ", sample_name)
  return(seurat_obj)
}

# -----------------------------
# EXECUTION LOOP
# -----------------------------
processed_objects <- list()

for(dir in cellranger_dirs){
  tryCatch({
    seurat_obj <- process_single_dataset(dir, species = species)
    if(!is.null(seurat_obj)) {
      processed_objects[[basename(dir)]] <- seurat_obj
    }
  }, error = function(e) {
    warning("Failed to process ", dir, ": ", e$message)
  })
}

# -----------------------------
# SUMMARY REPORT
# -----------------------------
if(length(processed_objects) > 0) {
  summary_data <- data.frame()
  for(sample_name in names(processed_objects)) {
    obj <- processed_objects[[sample_name]]
    summary_data <- rbind(summary_data, data.frame(
      Sample = sample_name,
      Cells = ncol(obj),
      Genes = nrow(obj),
      Clusters = length(unique(Idents(obj))),
      Mean_Features = mean(obj$nFeature_RNA),
      Mean_UMIs = mean(obj$nCount_RNA)
    ))
  }
  write.csv(summary_data, file.path(output_dir, "processing_summary.csv"), row.names = FALSE)
}

message("Pipeline execution completed. Processed ", length(processed_objects), " samples.")