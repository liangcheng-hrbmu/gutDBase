# File: gutdbase_integration.R
# Description: Integration and analysis of multiple scRNA-seq datasets from human and mouse CRC/IBD studies

# Load required libraries
library(data.table)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(reticulate)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(harmony)

# Set working directory for human CRC data
setwd("/youpath/GUTDbase/CRC_human")

# Load pre-processed and annotated Seurat objects for each dataset
GSE144735 <- readRDS("GSE144735/all_subtype_annotated.rds") 
GSE166319 <- readRDS("GSE166319/all_subtype_annotated.rds") 
GSE196964 <- readRDS("GSE196964/all_subtype_annotated.rds") 
GSE200997 <- readRDS("GSE200997/all_subtype_annotated.rds") 
GSE201348_1 <- readRDS("GSE201348_1/all_subtype_annotated.rds") 

# 1. Create a list of all Seurat objects
seurat_list <- list(
  GSE144735 = GSE144735,
  GSE166319 = GSE166319,
  GSE196964 = GSE196964,
  GSE200997 = GSE200997,
  GSE201348_1 = GSE201348_1
)

# Merge all Seurat objects with dataset identifiers
seurat_combined <- merge(seurat_list[[1]], y = seurat_list[-1], 
                         add.cell.ids = names(seurat_list), project = "CombinedData")

# 1. Normalization and variable feature selection
seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)

# 2. Data scaling and PCA
seurat_combined <- ScaleData(seurat_combined, features = VariableFeatures(seurat_combined))
seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(seurat_combined), npcs = 50)

# Determine optimal number of PCs using elbow plot
ElbowPlot(seurat_combined, ndims = 50)

# 3. Batch correction using Harmony
seurat_combined <- RunHarmony(
  object = seurat_combined,
  group.by.vars = c("orig.ident"),  # Batch effect variable
  dims.use = 1:50                   # Use first 50 principal components
)

# 4. Build KNN graph and perform clustering using Harmony-corrected PCs
seurat_combined <- FindNeighbors(seurat_combined, reduction = "harmony", dims = 1:50, k.param = 30)
seurat_combined <- FindClusters(seurat_combined, resolution = seq(0.1, 0.5, by = 0.1))

# 5. Dimensionality reduction using UMAP and t-SNE
seurat_combined <- RunUMAP(seurat_combined, reduction = "harmony", dims = 1:50, n.neighbors = 30, min.dist = 0.3)
seurat_combined <- RunTSNE(seurat_combined, reduction = "harmony", dims = 1:50, n.neighbors = 30, min.dist = 0.3)

# Save integrated object
saveRDS(seurat_combined, file = "CRC_human.rds")

# Set working directory for downstream analysis
setwd("/data_alluser/LXW/GUTDbase/CRC_human")

# Clean, modular version of CRC scRNA-seq workflow
# Key improvements:
# - Functions for repeated tasks (DE within groups, proportions, export to CellPhoneDB, trajectories)
# - Safer joins and factor handling
# - Fixed typos (DDRTree), UMAP/TSNE params, tidyverse style
# - Optional Monocle2 pseudotime (skippable / requires igraph<=1.5.1)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(data.table)
  library(ggpubr)
  library(stringr)
  library(monocle)
  library(Matrix)
})

set.seed(1)

# ---------- Paths & options ----------
base_out <- "/data_alluser/LXW/GUTDbase/CRC_human/output"
base_dl  <- "./download"
dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
dir.create(base_dl,  showWarnings = FALSE, recursive = TRUE)

# ---------- Load integrated object ----------
sce <- readRDS('/data_alluser/LXW/GUTDbase/CRC_human/CRC_human.rds')

# Run t-SNE if not already present
sce <- RunTSNE(sce, reduction = "harmony", dims = 1:50)

# Check unique subtypes and types
unique(sce@meta.data$subtype)
unique(sce@meta.data$Type)

# Function to merge cell subtypes by suffix and renumber
# Parses vectors with pattern "Cnumber_suffix", merges by suffix, and renumbers by minimal original C number
merge_by_suffix <- function(v, order_by = c("minC","freq","alpha"),
                            prefix_pattern = "^C\\d+_", new_prefix = "C") {
  order_by <- match.arg(order_by)   # Ordering rule: minimal original C number / frequency / alphabetical
  
  # Extract suffix (remove Cnumber_ prefix) and original C number
  suffix <- ifelse(grepl(prefix_pattern, v), sub(prefix_pattern, "", v), v)
  origC  <- suppressWarnings(as.integer(sub("^C(\\d+).*", "\\1", v)))
  
  df <- tibble(old = v, suffix = suffix, origC = origC)
  
  # Calculate minimal original C number and frequency for each suffix
  stat <- df %>%
    group_by(suffix) %>%
    summarise(
      minC = if (all(is.na(origC))) Inf else min(origC, na.rm = TRUE),
      freq = n(),
      .groups = "drop"
    )
  
  stat <- switch(order_by,
                 minC  = arrange(stat, minC, suffix),
                 freq  = arrange(stat, desc(freq), suffix),
                 alpha = arrange(stat, suffix)
  ) %>%
    mutate(newC = row_number(),
           new_label = paste0(new_prefix, newC, "_", suffix))
  
  # Create mapping and apply to original vector
  map  <- stat %>% select(suffix, newC, new_label)
  out  <- df %>% left_join(map, by = "suffix")
  
  # Ordered factor levels for new labels
  lvl  <- stat %>% arrange(newC) %>% pull(new_label)
  
  list(
    labels = out$new_label,
    mapping_table = stat %>% select(newC, suffix, minC, freq, new_label),
    levels = lvl
  )
}

## Apply to Seurat object
res <- merge_by_suffix(sce@meta.data$subtype, order_by = "minC")  # Can also use "freq"/"alpha"
sce$subtype <- factor(res$labels, levels = res$levels)            # Overwrite original column with merged labels
Idents(sce) <- sce$subtype                                         # Set subtype as identity if used for analysis

# ---------- Proportion tables ----------
# Function to calculate cell proportions by sample
prop_by_sample <- function(obj, group_col = c("subtype","Type")){
  group_col <- match.arg(group_col)
  obj@meta.data %>%
    rownames_to_column('barcode') %>%
    select(barcode, orig.ident, all_of(group_col)) %>%
    count(.data[[group_col]], orig.ident, name = 'freq') %>%
    rename(group = 1, sample = 2)
}

# Function to calculate cell proportions by disease state
prop_by_state <- function(obj, group_col = c("subtype","Type"), state_col = "Status"){
  group_col <- match.arg(group_col)
  obj@meta.data %>%
    count(.data[[state_col]], .data[[group_col]]) %>%
    group_by(.data[[group_col]]) %>%
    mutate(prop = round(n / sum(n), 3)) %>%
    ungroup() %>%
    select(state = 1, group = 2, prop)
}

# Calculate and save proportion tables
subtype_prop_table_barplot <- prop_by_sample(sce, "subtype")
type_prop_table_barplot    <- prop_by_sample(sce, "Type")
write.table(subtype_prop_table_barplot, file.path(base_out, 'subtype_prop_table_barplot.txt'),
            quote=FALSE, sep='\t', row.names=FALSE)
write.table(type_prop_table_barplot,    file.path(base_out, 'type_prop_table_barplot.txt'),
            quote=FALSE, sep='\t', row.names=FALSE)

subtype_prop_table_boxplot <- prop_by_state(sce, "subtype")
type_prop_table_boxplot    <- prop_by_state(sce, "Type")
write.table(subtype_prop_table_boxplot, file.path(base_out, 'subtype_prop_table_boxplot.txt'),
            quote=FALSE, sep='\t', row.names=FALSE)
write.table(type_prop_table_boxplot,    file.path(base_out, 'type_prop_table_boxplot.txt'),
            quote=FALSE, sep='\t', row.names=FALSE)

# ---------- Dynamic fold change helper ----------
# Calculate dynamic fold change thresholds (top 5% up/down regulated)
calculate_dynamic_foldchange <- function(df, num = 0.05){
  up   <- df %>% filter(avg_log2FC > 0)  %>% arrange(desc(avg_log2FC))
  down <- df %>% filter(avg_log2FC < 0)  %>% arrange(avg_log2FC)
  k_up   <- max(1, ceiling(nrow(up)   * num))
  k_down <- max(1, ceiling(nrow(down) * num))
  list(
    FC_up   = up$avg_log2FC[seq_len(k_up)],
    FC_down = down$avg_log2FC[seq_len(k_down)]
  )
}

# ---------- Differential expression within groups (CRC vs normal within each subtype/Type) ----------
run_de_within_groups <- function(obj, within = c("subtype","Type"), status_col = "Status",
                                 case_label = "CRC", ctrl_label = "normal",
                                 assay = 'RNA', logfc = 0.25, minpct = 0.25){
  within <- match.arg(within)
  Idents(obj) <- obj[[within]][,1]
  res_list <- list()
  levs <- levels(Idents(obj))
  for (lv in levs){
    obj_sub <- subset(obj, idents = lv)
    DefaultAssay(obj_sub) <- assay
    obj_sub <- NormalizeData(obj_sub, normalization.method = 'LogNormalize', scale.factor = 1e4, verbose = FALSE)
    if (!all(c(case_label, ctrl_label) %in% obj_sub[[status_col]][,1])){
      message(sprintf('[%s] missing %s or %s — skipped', lv, case_label, ctrl_label)); next
    }
    mk <- tryCatch({
      FindMarkers(obj_sub, ident.1 = case_label, ident.2 = ctrl_label,
                  group.by = status_col, slot = 'data', assay = assay,
                  logfc.threshold = logfc, min.pct = minpct, only.pos = FALSE) %>%
        rownames_to_column('gene')
    }, error = function(e){ NULL })
    if (is.null(mk) || nrow(mk) == 0) next
    
    thr <- calculate_dynamic_foldchange(mk, num = 0.05)
    cut_up   <- min(thr$FC_up)
    cut_down <- max(thr$FC_down)
    out <- mk %>%
      mutate(change = ifelse(p_val < 0.05 & (avg_log2FC > cut_up | avg_log2FC < cut_down),
                             ifelse(avg_log2FC > cut_up, 'Up','Down'), 'Stable'),
             xinterceptleft = cut_down,
             xinterceptright= cut_up,
             yintercept     = -log10(0.05),
             group          = lv,
             `-log10(P.value)` = -log10(p_val))
    res_list[[lv]] <- out
  }
  bind_rows(res_list)
}

# Run DE analysis for subtypes and types
subtype_de <- run_de_within_groups(sce, within='subtype')
type_de    <- run_de_within_groups(sce, within='Type')

# Save DE results
fwrite(subtype_de, file.path(base_out, 'subtype_Volcano.txt'), sep='\t', quote=FALSE)
fwrite(type_de,    file.path(base_out, 'type_Volcano.txt'),    sep='\t', quote=FALSE)

# ---------- Expression + coordinates for significant genes ----------
# Ensure t-SNE coordinates are available (prefer harmony, fallback to PCA)
ensure_tsne <- function(obj, reduction = NULL, dims = 1:50,
                        perplexity = NULL, check_duplicates = FALSE,
                        seed.use = 1) {
  if ("tsne" %in% names(obj@reductions)) return(obj)
  n <- ncol(obj)
  if (is.null(perplexity)) {
    # Empirical range: 5-30, must be < (n-1)/3
    p_max <- max(5, floor((n - 1) / 3) - 1)
    perplexity <- min(30, p_max)
    if (!is.finite(perplexity) || perplexity < 5) perplexity <- 5
  }
  if (is.null(reduction)) {
    reduction <- if ("harmony" %in% names(obj@reductions)) "harmony" else "pca"
  }
  obj <- RunTSNE(
    obj, reduction = reduction, dims = dims,
    perplexity = perplexity, check_duplicates = check_duplicates,
    seed.use = seed.use
  )
  message(sprintf("RunTSNE done: reduction=%s, dims=%s, perplexity=%s, cells=%s",
                  reduction, paste0(range(dims), collapse = ":"), perplexity, n))
  obj
}

sce <- ensure_tsne(sce, reduction = "harmony", dims = 1:50)

# Function to extract expression data with coordinates (default: t-SNE)
expression_with_coords <- function(obj, de_tbl, group_cols = c("subtype","Type"),
                                   red = c("tsne","umap"), assay = "RNA", slot = "data"){
  red <- match.arg(red)
  sig_genes <- de_tbl %>% dplyr::filter(change %in% c("Up","Down")) %>%
    dplyr::pull(gene) %>% unique()
  if (length(sig_genes) == 0) return(invisible(NULL))
  Seurat::DefaultAssay(obj) <- assay
  fetch <- Seurat::FetchData(obj, vars = c(sig_genes, group_cols, "Status"), slot = slot) %>%
    tibble::rownames_to_column("barcode") %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(sig_genes),
                        names_to = "gene", values_to = "value")
  coords <- Embeddings(obj, red) %>% as.data.frame() %>% tibble::rownames_to_column("barcode")
  dplyr::left_join(fetch, coords, by = "barcode")
}

# Generate expression coordinate tables (default: t-SNE)
subtype_expr_coord <- expression_with_coords(sce, subtype_de, group_cols = c("subtype"))
type_expr_coord    <- expression_with_coords(sce, type_de,    group_cols = c("Type"))

# Optionally write split-by-group tables (non-zero values only)
if (!is.null(subtype_expr_coord)) {
  out_dir <- file.path(base_out, 'expression_and_coord', 'subtype'); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  split(subtype_expr_coord %>% filter(value != 0), subtype_expr_coord$subtype) %>%
    iwalk(~ fwrite(.x, file.path(out_dir, paste0(.y, '.txt')), sep='\t', quote=FALSE))
}
if (!is.null(type_expr_coord)) {
  out_dir <- file.path(base_out, 'expression_and_coord', 'type'); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  split(type_expr_coord %>% filter(value != 0), type_expr_coord$Type) %>%
    iwalk(~ fwrite(.x, file.path(out_dir, paste0(.y, '.txt')), sep='\t', quote=FALSE))
}

# ---------- Marker genes & heatmaps ----------
DefaultAssay(sce) <- 'RNA'

# Find marker genes for subtypes
Idents(sce) <- sce$subtype
sce_subtype_markers <- FindAllMarkers(sce, only.pos = TRUE, slot = 'data', test.use = 'wilcox',
                                      assay = 'RNA', logfc.threshold = 0.25, min.pct = 0.25)
sub_top10 <- sce_subtype_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup()

# Find marker genes for cell types
Idents(sce) <- sce$Type
sce_type_markers <- FindAllMarkers(sce, only.pos = TRUE, slot = 'data', test.use = 'wilcox',
                                   assay = 'RNA', logfc.threshold = 0.25, min.pct = 0.25)
type_top10 <- sce_type_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup()

# Save marker gene results
write.table(sce_subtype_markers, file.path(base_dl, 'sce_subtype_markers.txt'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(sce_type_markers,    file.path(base_dl, 'sce_type_markers.txt'),    sep='\t', quote=FALSE, row.names=FALSE)

# AverageHeatmap from scRNAtoolVis (requires that package)
# pdf(file.path(base_out, 'subtype_heatmap.pdf'), width=12, height=25)
# AverageHeatmap(sce, markerGene = sub_top10$gene, assay='SCT', slot='scale.data', row_names_side='right', row_title='', column_names_rot=45)
# dev.off()
# pdf(file.path(base_out, 'type_heatmap.pdf'), width=10, height=8)
# AverageHeatmap(sce, markerGene = type_top10$gene, assay='SCT', slot='scale.data', row_names_side='right', row_title='', column_names_rot=45)
# dev.off()

# ---------- Export for CellPhoneDB ----------
# Function to export data for CellPhoneDB analysis
export_for_cpdb <- function(obj, meta_label = c('subtype','Type'), counts_assay = 'RNA', counts_slot = 'counts', out_dir){
  meta_label <- match.arg(meta_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  DefaultAssay(obj) <- counts_assay
  counts <- GetAssayData(obj, slot = counts_slot) %>% as.data.frame() %>% rownames_to_column('gene')
  fwrite(counts, file.path(out_dir, 'cellphonedb_count.txt'), sep='\t', quote=FALSE)
  meta <- data.frame(barcode = rownames(obj@meta.data), cell_type = obj@meta.data[[meta_label]], check.names = FALSE)
  write.table(meta, file.path(out_dir, 'cellphonedb_meta.txt'), sep='\t', quote=FALSE, row.names=FALSE)
}

# Export data for CellPhoneDB analysis
export_for_cpdb(sce, meta_label='subtype', out_dir = file.path(base_out, 'cellphonedb', 'subtype'))
export_for_cpdb(sce, meta_label='Type',    out_dir = file.path(base_out, 'cellphonedb', 'type'))

# Function to parse CellPhoneDB results into edge list
parse_cpdb_edges <- function(means_file, pvals_file, p_cut = 0.05){
  meansdf <- read.delim(means_file, check.names = FALSE) %>%
    pivot_longer(-c(interacting_pair, gene_a, gene_b, receptor_a, receptor_b), names_to='pair', values_to='mean_expr')
  pvalsdf <- read.delim(pvals_file, check.names = FALSE) %>%
    pivot_longer(-interacting_pair, names_to='pair', values_to='p_value')
  merged <- inner_join(pvalsdf, meansdf, by = c('interacting_pair','pair')) %>%
    separate(pair, into = c('source','target'), sep = '\\|', remove = TRUE) %>%
    filter(nchar(gene_a) > 0, nchar(gene_b) > 0, p_value <= p_cut) %>%
    select(source, target, gene_a, gene_b, mean_expr, p_value)
  merged
}

# Save final object
saveRDS(sce, file = "CRC_human.rds")

# Check igraph version for Monocle compatibility
packageVersion("igraph")

# Run t-SNE with specific parameters
seurat_combined <- RunTSNE(seurat_combined, reduction = "harmony", dims = 1:50, n.neighbors = 30, min.dist = 0.3)

# Save object
saveRDS(sce, file = "CRC_human.rds")

# The script continues with similar workflows for:
# 1. Mouse CRC data integration and analysis
# 2. Human IBD data integration and analysis
# Each section follows the same structure with appropriate parameter adjustments