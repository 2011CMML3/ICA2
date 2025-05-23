if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("https://github.com/MarcElosua/SPOTlight")

library(SPOTlight)
library(Seurat)      
library(NMF)         
library(ggplot2)    

# --------------------------
# 1. Load Data
# --------------------------
# Load single-cell reference data (Seurat object)
# Should contain:
# - Normalized counts in RNA assay
# - Cell type labels in metadata column 'cell_type'
sc_data <- readRDS("sc_data.rds")

# Load spatial transcriptomics data (Seurat object)
# Should contain:
# - Normalized counts in RNA assay
# - Spatial coordinates in @images slot
spatial_data <- readRDS("spatial_data.rds")

# --------------------------
# 2. Preprocess Data
# --------------------------
# Normalize and find variable features for single-cell data
sc_data <- NormalizeData(sc_data)
sc_data <- FindVariableFeatures(sc_data)

# Normalize spatial data (using same features)
spatial_data <- NormalizeData(spatial_data)
spatial_data <- ScaleData(spatial_data, features = VariableFeatures(sc_data))

# Ensure common genes
common_genes <- intersect(rownames(sc_data), rownames(spatial_data))
sc_data <- subset(sc_data, features = common_genes)
spatial_data <- subset(spatial_data, features = common_genes)

# --------------------------
# 3. Run SPOTlight Deconvolution
# --------------------------
nmf_model <- SPOTlight(
  x = GetAssayData(sc_data, slot = "data"),  # Normalized scRNA-seq data
  y = GetAssayData(spatial_data, slot = "data"),  # Normalized spatial data
  groups = sc_data$cell_type,  # Cell type labels
  hvg = VariableFeatures(sc_data),  # Highly variable genes
  K = 10,  # Number of factors (default: auto-detected)
  min_cont = 0.01  # Minimum contribution threshold
)

# --------------------------
# 4. Extract and Analyze Results
# --------------------------
# Get cell type proportions matrix (spots x cell types)
cell_proportions <- nmf_model$weights
head(cell_proportions)

# Add results to spatial data metadata
spatial_data$spotlight_deconv <- cell_proportions

# Get NMF model details
nmf_details <- nmf_model$NMF

# --------------------------
# 5. Save Results
# --------------------------
# Save cell type proportions
write.csv(cell_proportions, "spotlight_cell_proportions.csv")

# Save full model
saveRDS(nmf_model, "spotlight_model.rds")

# Save spatial data with deconvolution results
saveRDS(spatial_data, "spatial_data_with_spotlight.rds")
