# bulk-human/bulk-mouse
# Using human IBS disease as example

#### Step 1: Batch correction by GSE number ####
library(GEOquery)
library(dplyr)
library(data.table)
library(limma)
library(gmodels)
library(ggpubr)
library(ggthemes)
library(tidyverse)
library(org.Hs.eg.db)
library(simpleaffy)
library(RColorBrewer)
library(magrittr)
library(DESeq2)
library(org.Mm.eg.db) 

# Set working directory and load data
setwd('/youpath/SingleCellData/Gut_Database/IBS/human/total_RNA/counts')
counts1 <- read.table("/youpath/SingleCellData/Gut_Database/IBS/human/total_RNA/counts/merged_mat.txt", 
                      header = TRUE, row.names = NULL, sep = "\t")
p_data <- read.table("/youpath/SingleCellData/Gut_Database/IBS/human/total_RNA/counts/p_data.txt", 
                     header = TRUE, row.names = NULL, sep = "\t")

# Add accession information to p_data
library(sva)
p_data$accession <- c(
  rep("GSE13367", 38),          # Rows 1-38
  rep("GSE14841", 9),           # Rows 39-47
  rep("GSE36701", 127),         # Rows 48-174
  rep("GSE63379", 67)           # Rows 175-241
)

# Save updated p_data
p_data %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/p_data.txt',
              row.names = T, sep = '\t', quote = F)

# For microarray data use this approach
counts1 <- column_to_rownames(counts1, "row.names")

# Build batch information (by GSE_id)
batch <- p_data$accession  # Batch = GSE13367/GSE14841/GSE36701/GSE63379

# Build group model (case/control, preserving biological differences)
mod <- model.matrix(~ factor(metaclass), data = p_data)

# Batch correction (ensure expression matrix column names match p_data$sample_id)
combat_data <- ComBat(
  dat = counts1,
  batch = batch,
  mod = mod
)

# Save batch-corrected data
combat_data %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/combat_data.txt',
              row.names = T, sep = '\t', quote = F)

#################
# For high-throughput sequencing data use this approach
###############

# Load required packages
library(sva)
library(dplyr)
library(tibble)

# 1. Process expression matrix: preserve sample IDs and confirm order
expr_sample_df <- tibble(
  sample = colnames(counts1),     # Extract sample IDs from expression matrix (fixed order)
  expr_order = 1:ncol(counts1)    # Mark sample order in expression matrix (for alignment)
)
cat("Expression matrix sample order (first 5):\n")
print(head(expr_sample_df))

# 2. Process clinical data p_data: bind sample IDs, align with expression matrix order
p_data_aligned <- p_data %>%
  dplyr::select(sample, accession, metaclass) %>%  # Keep only needed columns
  inner_join(expr_sample_df, by = "sample") %>%   # Keep only samples present in expression matrix
  arrange(expr_order) %>%  # Force sort by expression matrix sample order
  dplyr::select(-expr_order)  # Remove order marker column

# 3. Extract aligned batch and group information (order matches expression matrix exactly)
batch <- p_data_aligned$accession  # Batch information (order = expression matrix sample order)
metaclass_aligned <- p_data_aligned$metaclass  # Group information (same order)

# 4. Build group model mod: based on aligned group information, row order = sample order
mod <- model.matrix(
  ~ factor(metaclass_aligned),  # Use aligned group information
  data = data.frame(metaclass_aligned = metaclass_aligned)  # Data frame row order = sample order
)
# Add row names to mod (sample IDs), ensure correspondence with batch samples
rownames(mod) <- p_data_aligned$sample

# 5. Final validation: sample count, order, row names must be identical
cat("\n=== Final validation ===\n")
cat("\nExpression matrix sample count: ", ncol(counts1))
cat("\nBatch sample count: ", length(batch))
cat("\nMod row count: ", nrow(mod))
cat("\nBatch sample order (first 5): ", paste(head(batch), collapse = ", "))
cat("\nMod row names (first 5, sample IDs): ", paste(head(rownames(mod)), collapse = ", "))
cat("\nExpression matrix column names (first 5, sample IDs): ", paste(head(colnames(counts1)), collapse = ", "))

# If any inconsistencies exist, report error and show differences
if (!all(rownames(mod) == colnames(counts1)) || !all(p_data_aligned$sample == colnames(counts1))) {
  stop("Sample order/ID mismatch!\n",
       "Mod row names: ", paste(head(rownames(mod)), "..."), "\n",
       "Expression matrix column names: ", paste(head(colnames(counts1)), "..."))
}

# Additional validation and processing for ComBat_seq
cat("counts1 original type: ", class(counts1), "\n")
# If data.frame, convert to matrix; if non-numeric, force conversion
counts1 <- as.matrix(counts1)
# Check for non-numeric values (NA/characters) and handle
if (any(is.na(counts1))) {
  counts1[is.na(counts1)] <- 0  # Replace NA values with 0 (common for RNA-seq counts)
  cat("Replaced NA values in counts1 with 0\n")
}
if (!is.numeric(counts1)) {
  counts1 <- apply(counts1, 2, as.numeric)  # Force conversion to numeric
  cat("Forced conversion of counts1 to numeric matrix\n")
}
cat("Processed counts1 type: ", class(counts1), "\n")
cat("counts1 first 3 rows and 3 columns preview:\n")
print(head(counts1[, 1:3], 3))

# Simplify group model mod (avoid collinearity issues)
cat("\n=== Rebuilding group model ===\n")
# Option 1: Use group vector instead of design matrix, simpler and less error-prone
group_vector <- as.factor(p_data_aligned$metaclass)  # Directly use group vector (case/control)
cat("\nGroup vector levels (group types): ", levels(group_vector), "\n")
cat("Group vector length: ", length(group_vector), "\n")

# Final consistency validation
cat("\n=== Final consistency validation ===\n")
cat("\ncounts1 column count (samples): ", ncol(counts1))
cat("\nBatch length (samples): ", length(batch))
cat("\nGroup_vector length (samples): ", length(group_vector))
cat("\ncounts1 column names first 5: ", paste(head(colnames(counts1)), collapse = ", "))
cat("\np_data_aligned sample names first 5: ", paste(head(p_data_aligned$sample), collapse = ", "))
cat("\ngroup_vector first 5 groups: ", paste(head(group_vector), collapse = ", "))

# Confirm no sample order misalignment
if (!all(colnames(counts1) == p_data_aligned$sample)) {
  stop("Expression matrix column names don't match clinical data sample names!")
}

# Execute ComBat_seq (using group_vector)
cat("\n=== Starting ComBat-seq ===\n")
combat_data <- ComBat_seq(
  counts = counts1,        # Numeric matrix (verified)
  batch = batch,           # Batch vector (aligned)
  group = group_vector     # Group vector (more stable than design matrix)
)

# Save results
write.table(
  combat_data,
  'E:/Desktop/combat_data.txt',
  row.names = TRUE,
  sep = '\t',
  quote = FALSE
)

cat("\n=== Batch correction completed ===\n")
cat("\nPost-batch correction matrix dimensions: ", dim(combat_data), "\n")
cat("Post-batch correction matrix first 3 rows and 3 columns preview:\n")
print(head(combat_data[, 1:3], 3))

#################
#### Split total expression matrix by accession ####
library(dplyr)
library(tibble)

# 1. Extract "sample ID-accession" mapping from p_data (deduplicate)
sample_accession_map <- p_data %>%
  dplyr::select(sample, accession) %>%  # Keep only needed columns
  distinct(sample, .keep_all = TRUE) %>%  # Deduplicate
  filter(sample %in% colnames(combat_data))  # Keep only samples in combat_data

# 2. Check accession distribution
table(sample_accession_map$accession)

# Extract unique accession values
unique_accessions <- unique(sample_accession_map$accession)

# 3. Split combat_data by accession (generate data frames stored in list)
combat_split_list <- lapply(unique_accessions, function(acc) {
  # Step A: Filter sample IDs for current accession
  samples_for_acc <- sample_accession_map %>%
    dplyr::filter(accession == acc) %>%
    dplyr::pull(sample)  # Extract sample ID vector
  
  # Step B: Filter these sample columns from combat_data
  combat_subset <- combat_data %>%
    as.data.frame() %>%  # Ensure data frame format
    dplyr::select(all_of(samples_for_acc))  # Filter columns by sample ID
  
  # Step C: Add gene name column
  combat_subset <- combat_subset %>%
    tibble::rownames_to_column(var = "gene")
  
  return(combat_subset)
})

# 4. Name the data frames in the list
names(combat_split_list) <- unique_accessions

# Extract individual datasets
combat_1 <- combat_split_list[["GSE13367"]]  # First accession
combat_2 <- combat_split_list[["GSE14841"]]  # Second accession
combat_3 <- combat_split_list[["GSE36701"]]  # Third accession
combat_4 <- combat_split_list[["GSE63379"]]  # Fourth accession

######## Step 2: Differential expression analysis ####
p_data <- read.table("/youpath/SingleCellData/Gut_Database/IBS/human/total_RNA/GSE63379/download/p_data.txt", 
                     header = TRUE, row.names = NULL, sep = "\t")
combat_4 <- column_to_rownames(combat_4, "gene")
exp <- combat_4
metaclass <- factor(p_data$metaclass) # Clinical data

# Design matrix for linear modeling
design <- model.matrix(~0 + factor(metaclass))
colnames(design) <- levels(metaclass)
rownames(design) <- colnames(exp)

# Create contrast matrix
contrast.matrix <- makeContrasts(paste0(rev(unique(metaclass)), collapse = '-'), levels = design) # control-case

# Fit linear model
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 

# Extract results
tempOutput = topTable(fit2, coef = 1, n = Inf)
nrDEG = na.omit(tempOutput) %>% rownames_to_column(., var = 'gene')

# Function to calculate dynamic fold change thresholds
calculate_dynamic_foldchange <- function(df, num){
  up_gene_diff <- length(which(df$logFC > 0))
  up_gene_pct <- ceiling(up_gene_diff * num)
  
  down_gene_diff <- length(which(df$logFC < 0))
  down_gene_pct <- ceiling(down_gene_diff * num)
  
  gene_up <- df[which(df$logFC > 0), ]
  gene_up_pct <- gene_up[order(-gene_up$logFC), ][1:up_gene_pct, 'gene']
  
  gene_down <- df[which(df$logFC < 0), ]
  gene_down_pct <- gene_down[order(gene_down$logFC), ][1:down_gene_pct, 'gene']
  
  All_gene_pct <- df[which(df$gene %in% c(gene_up_pct, gene_down_pct)), ]
  
  FC_up <- All_gene_pct[which(All_gene_pct$logFC > 0), ]$logFC
  FC_down <- All_gene_pct[which(All_gene_pct$logFC < 0), ]$logFC
  
  return(list(FC_up, FC_down))
}

# Calculate cutoff values
cut_off_logFC_up <- min(calculate_dynamic_foldchange(nrDEG, num = 0.05)[[1]])   # Dynamic fold change threshold
cut_off_logFC_down <- max(calculate_dynamic_foldchange(nrDEG, num = 0.05)[[2]])

# Classify genes as Up, Down, or Stable
nrDEG %>% 
  mutate(group = case_when(
    logFC >= cut_off_logFC_up & P.Value <= 0.05 ~ "Up",
    logFC <= cut_off_logFC_down & P.Value <= 0.05 ~ "Down",
    TRUE ~ "Stable"
  )) -> res_4

# Filter and save significant results
res_4_filtered <- res_4[res_4$group %in% c("Up", "Down"), ]
res_4_filtered %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/res_4_deg.txt',
              row.names = T, sep = '\t', quote = F)

######### Volcano plot data preparation
xinterceptleft = cut_off_logFC_down
xinterceptright = cut_off_logFC_up
yintercept = -log10(0.05)

volcano1 <- data.frame(yintercept = yintercept,
                       xinterceptleft = rep(xinterceptleft, nrow(res_4)),
                       xinterceptright = rep(xinterceptright, nrow(res_4)),
                       p.val = -log10(res_4$P.Value),
                       res_4) %>%
  dplyr::select(-c('AveExpr', 't', 'B', 'adj.P.Val')) %>%
  dplyr::rename(`-log10(p.val)` = p.val)

volcano1 %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/hum_volcanodata_4.txt',
              row.names = T, sep = '\t', quote = F)

######### Step 3: Enrichment analysis ####
####### Split and analyze up/down regulated genes separately
########## GO enrichment analysis
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db) 

# Extract up and down regulated genes
up_genes <- res_4 %>% 
  filter(group == "Up") %>% 
  pull(gene) %>% 
  unique()

down_genes <- res_4 %>% 
  filter(group == "Down") %>% 
  pull(gene) %>% 
  unique()

# Check differential gene counts
cat("Up-regulated gene count:", length(up_genes), "\n")
cat("Down-regulated gene count:", length(down_genes), "\n")

# Gene ID conversion (Symbol -> EntrezID)
# Up-regulated genes
up_entrez <- bitr(
  geneID = up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db # For mouse: org.Mm.eg.db
) %>% pull(ENTREZID) %>% unique()

# Down-regulated genes
down_entrez <- bitr(
  geneID = down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db # For mouse: org.Mm.eg.db
) %>% pull(ENTREZID) %>% unique()

# Check converted counts
cat("Converted up-regulated gene ID count:", length(up_entrez), "\n")
cat("Converted down-regulated gene ID count:", length(down_entrez), "\n")

# Save gene lists
up_entrez %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/up4.txt',
              row.names = T, sep = '\t', quote = F)

down_entrez %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/down4.txt',
              row.names = T, sep = '\t', quote = F)

# 2.1 Up-regulated genes GO enrichment
go_up <- clusterProfiler::enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db, # For mouse: org.Mm.eg.db
  keyType = "ENTREZID",
  ont = "ALL",  # Analyze BP/CC/MF simultaneously
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# 2.2 Down-regulated genes GO enrichment
go_down <- clusterProfiler::enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db, # For mouse: org.Mm.eg.db
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# 2.3 Organize GO results (add status column and standardize format)
go_up_df <- as.data.frame(go_up) %>% 
  mutate(
    status = "Up",
    logpvalue = -log10(pvalue)
  ) %>% 
  dplyr::select(
    ontology = ONTOLOGY,
    goid = ID,
    description = Description,
    pvalue,
    logpvalue,
    geneid = geneID,
    status
  )

go_down_df <- as.data.frame(go_down) %>% 
  mutate(
    status = "Down",
    logpvalue = -log10(pvalue)
  ) %>% 
  dplyr::select(
    ontology = ONTOLOGY,
    goid = ID,
    description = Description,
    pvalue,
    logpvalue,
    geneid = geneID,
    status
  )

# Combine up and down GO results
hum_go <- dplyr::bind_rows(go_up_df, go_down_df)
hum_go %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/hum_godown4.txt',
              row.names = T, sep = '\t', quote = F)

######## KEGG enrichment analysis (local)
library(clusterProfiler)
up_entrez <- as.character(read.table("E:/Desktop/up1.txt", header = TRUE, row.names = NULL, sep = "\t")[, 2])
down_entrez <- as.character(read.table("E:/Desktop/down1.txt", header = TRUE, row.names = NULL, sep = "\t")[, 2])

# 3.1 Up-regulated genes KEGG enrichment
kegg_up <- enrichKEGG(
  gene = up_entrez,
  organism = "hsa",  # Human: hsa, Mouse: mmu
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2
)

# 3.2 Down-regulated genes KEGG enrichment
kegg_down <- enrichKEGG(
  gene = down_entrez,
  organism = "hsa", # Human: hsa, Mouse: mmu
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2
)

######### Up-regulated
# 1. Read EntrezID and create data frame
up_entrez_df <- data.frame(
  entrez_id = as.character(read.table(
    "E:/Desktop/up1.txt", 
    header = TRUE, 
    row.names = NULL,
    sep = "\t"
  )[, 2])  # Extract second column as EntrezID
)

# 2. Add Symbol column (via org.Hs.egSYMBOL conversion)
up_entrez_df$symbol <- as.character(
  mget(up_entrez_df$entrez_id, org.Hs.egSYMBOL, ifnotfound = NA)
) # For mouse: org.Mm.egSYMBOL

# 3.3 Organize KEGG results (add status column and standardize format)
kegg_up_df <- as.data.frame(kegg_up) %>% 
  mutate(
    ontology = "KEGG",
    status = "Up",
    logpvalue = -log10(pvalue),
    # Split geneID and match Symbol
    gene_symbol = sapply(strsplit(as.character(geneID), "/"), function(ids) {
      # Use data frame matching: find Symbol for each EntrezID
      symbols <- up_entrez_df$symbol[match(ids, up_entrez_df$entrez_id)]
      paste(ifelse(is.na(symbols), "NA", symbols), collapse = "/")
    })
  ) %>% 
  dplyr::select(
    ontology,
    goid = ID,
    description = Description,
    pvalue,
    logpvalue,
    geneid = gene_symbol,
    status
  )

######### Down-regulated
down_entrez_df <- data.frame(
  entrez_id = as.character(read.table(
    "E:/Desktop/down1.txt", 
    header = TRUE, 
    row.names = NULL,
    sep = "\t"
  )[, 2])  # Extract second column as EntrezID
)

# 2. Add Symbol column (via org.Hs.egSYMBOL conversion)
down_entrez_df$symbol <- as.character(
  mget(down_entrez_df$entrez_id, org.Hs.egSYMBOL, ifnotfound = NA)
) # For mouse: org.Mm.egSYMBOL

# 3.3 Organize KEGG results (add status column and standardize format)
kegg_down_df <- as.data.frame(kegg_down) %>% 
  mutate(
    ontology = "KEGG",
    status = "down",
    logpvalue = -log10(pvalue),
    # Split geneID and match Symbol
    gene_symbol = sapply(strsplit(as.character(geneID), "/"), function(ids) {
      # Use data frame matching: find Symbol for each EntrezID
      symbols <- down_entrez_df$symbol[match(ids, down_entrez_df$entrez_id)]
      paste(ifelse(is.na(symbols), "NA", symbols), collapse = "/")
    })
  ) %>% 
  dplyr::select(
    ontology,
    goid = ID,
    description = Description,
    pvalue,
    logpvalue,
    geneid = gene_symbol,
    status
  )

# Combine up and down KEGG results
hum_kegg <- bind_rows(kegg_up_df, kegg_down_df)
hum_kegg %>%
  write.table(., 'E:/Desktop/hum_keggdown1.txt',
              row.names = T, sep = '\t', quote = F)

######### Step 4: PCA analysis ####
p_data <- read.table("/youpath/SingleCellData/Gut_Database/IBS/human/total_RNA/counts/p_data.txt", 
                     header = TRUE, row.names = NULL, sep = "\t")
exp <- combat_data
metaclass <- factor(p_data$metaclass) # Clinical data

# Design matrix for linear modeling
design <- model.matrix(~0 + factor(metaclass))
colnames(design) <- levels(metaclass)
rownames(design) <- colnames(exp)

# Create contrast matrix
contrast.matrix <- makeContrasts(paste0(rev(unique(metaclass)), collapse = '-'), levels = design) # control-case

# Fit linear model
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 

# Extract results
tempOutput = topTable(fit2, coef = 1, n = Inf)
nrDEG = na.omit(tempOutput) %>% rownames_to_column(., var = 'gene')

# Calculate dynamic fold change thresholds
cut_off_logFC_up <- min(calculate_dynamic_foldchange(nrDEG, num = 0.05)[[1]])   # Dynamic fold change threshold
cut_off_logFC_down <- max(calculate_dynamic_foldchange(nrDEG, num = 0.05)[[2]])

# Classify genes
nrDEG %>% 
  mutate(group = case_when(
    logFC >= cut_off_logFC_up & P.Value <= 0.05 ~ "Up",
    logFC <= cut_off_logFC_down & P.Value <= 0.05 ~ "Down",
    TRUE ~ "Stable"
  )) -> res

# Function to calculate PCA with expression data
calculate_pca_with_expr <- function(res_2, counts, pdata) {
  # Step 1: Extract expression data (exclude gene column) + gene names
  expr_data <- counts %>% dplyr::select(-gene)  # Sample expression values (numeric)
  gene_names <- counts$gene             # Gene name list
  
  # Step 2: Calculate variance for each row (gene) using base R
  gene_variance <- apply(
    X = expr_data,
    MARGIN = 1,  # 1 = calculate by row (each row = 1 gene)
    FUN = function(row) var(row, na.rm = TRUE)
  )
  
  # Step 3: Sort by variance descending, select top 500 genes
  variance_df <- data.frame(
    gene = gene_names,
    var = gene_variance,
    stringsAsFactors = FALSE
  ) %>% dplyr::arrange(dplyr::desc(var))  # Sort variance from high to low
  selected_gene_num <- min(500, nrow(variance_df))  # Avoid insufficient genes
  selected_genes <- variance_df$gene[1:selected_gene_num]
  
  # Step 4: Extract expression data for selected genes, transpose
  selected_expr <- expr_data[counts$gene %in% selected_genes, , drop = FALSE]
  pca_input <- t(selected_expr)  # Transpose: rows = samples, columns = genes
  
  # Step 5: Perform PCA
  pca <- prcomp(pca_input, scale. = FALSE)
  
  # Step 6: Organize PCA results (coordinates + variance + clinical info)
  pca_coords <- as.data.frame(pca$x[, 1:2])  # PC1 and PC2 coordinates
  colnames(pca_coords) <- c("PC1", "PC2")
  pca_coords$sample <- rownames(pca_coords)  # Sample ID
  
  # Calculate principal component contribution (variance explained %)
  total_var <- sum(pca$sdev^2)
  pc1_pct <- round((pca$sdev[1]^2 / total_var) * 100, 2)
  pc2_pct <- round((pca$sdev[2]^2 / total_var) * 100, 2)
  
  # Merge clinical data (metaclass grouping)
  pca_data <- merge(
    x = pca_coords,
    y = p_data[, c("sample", "metaclass")],
    by = "sample",
    all.x = TRUE
  ) %>% dplyr::mutate(
    PC1_var = paste0(pc1_pct, "%"),
    PC2_var = paste0(pc2_pct, "%")
  )
  
  # Step 7: Merge differential gene expression data
  degs <- res_2 %>% dplyr::filter(group %in% c("Up", "Down"))
  expr_to_merge <- if (nrow(degs) == 0) {
    counts %>% dplyr::filter(gene %in% selected_genes)  # Use top 500 genes if no DEGs
  } else {
    counts %>% dplyr::filter(gene %in% degs$gene)       # Use DEGs if available
  }
  
  # Convert to long format + merge PCA information
  final_dt <- expr_to_merge %>%
    tidyr::pivot_longer(cols = !gene, names_to = "sample", values_to = "expression") %>%
    merge(x = ., y = pca_data, by = "sample", all = FALSE)
  
  return(final_dt)
}

combat_data <- as.data.frame(combat_data)
combat_data <- rownames_to_column(combat_data, var = "gene")

# Calculate PCA coordinates
exp_pca_Degs_coord <- calculate_pca_with_expr(res_2 = res,
                                              counts = combat_data,
                                              pdata = p_data)

exp_pca_Degs_coord %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/hum_pca.txt',
              row.names = T, sep = '\t', quote = F)

######### Step 5: Immune infiltration analysis ####
library(tidyverse)
library(rstatix)
library(ggpubr)
library(CIBERSORT)

combat_data <- read.table("/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/combat_data.txt", 
                           header = TRUE, row.names = NULL, sep = "\t")
p_data <- read.table("/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/p_data.txt", 
                     header = TRUE, row.names = NULL, sep = "\t")

# Prepare data for CIBERSORT
gene <- rownames(combat_data)
train_dat <- cbind(gene, combat_data)

# Ensure first column name is "Gene"
train_dat <- train_dat[, -1]
colnames(train_dat)[1] <- "Gene"
train_dat %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/train_dat.txt',
              row.names = T, sep = '\t', quote = F)

# Run CIBERSORT
source("/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/CIBERSORT.R")
LM22.txt <- "/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/LM22.txt"

results <- CIBERSORT('/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/LM22.txt',
                     '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/train_dat.txt', 
                     perm = 1000, QN = T)

# Process results
results <- as.data.frame(results)
results <- dplyr::select(results, -c('P-value', 'Correlation', 'RMSE'))
results$sample <- rownames(results)

# Convert to long format
longer_data <- results %>% 
  pivot_longer(
    cols = -c("sample"),
    names_to = "ImmuneCell",
    values_to = "Score"
  )

# Merge with clinical data
longer_data <- merge(longer_data, p_data, by = "sample")
longer_data %>%
  write.table(., '/youpath/SingleCellData/Gut_Database/batch_correction/human/IBS/Expression_profiling_by_array/hum_box.txt',
              row.names = T, sep = '\t', quote = F)