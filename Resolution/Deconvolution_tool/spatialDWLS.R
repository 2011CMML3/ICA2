if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("RubD/Giotto")

# Load required libraries
library(Giotto)
library(Seurat)  # For data handling (optional)

# --------------------------
# 1. Load Data
# --------------------------
# Load single-cell reference data (Seurat object format)
sc_data <- readRDS("sc_data.rds") 

# Load spatial transcriptomics data (Seurat object format)
spatial_data <- readRDS("spatial_data.rds")  

# --------------------------
# 2. Preprocess Data
# --------------------------
# Convert to Giotto objects if not already
# giotto_sc <- createGiottoObject(raw_exprs = sc_data@assays$RNA@counts)
# giotto_spatial <- createGiottoObject(raw_exprs = spatial_data@assays$RNA@counts)

# Ensure both datasets use the same genes
common_genes <- intersect(rownames(sc_data), rownames(spatial_data))
sc_data <- sc_data[common_genes, ]
spatial_data <- spatial_data[common_genes, ]

# Normalize data (log-normalization recommended)
sc_data <- NormalizeData(sc_data)
spatial_data <- NormalizeData(spatial_data)

# --------------------------
# 3. Run spatialDWLS Deconvolution
# --------------------------
dwls_results <- runDWLSDeconv(
  spat_obj = spatial_data,      
  sc_obj = sc_data,              
  cell_type_column = "cell_type", 
  sig_gene_method = "DWLS",      
  nr_permutations = 100,         
  return_spatial_obj = FALSE      
)

# --------------------------
# 4. Extract and Analyze Results
# --------------------------
# Get cell type proportions matrix (spots x cell types)
cell_proportions <- dwls_results$weights
head(cell_proportions)

# Add results to spatial data metadata
spatial_data$dwls_deconv <- cell_proportions



# --------------------------
# 6. Save Results
# --------------------------
# Save deconvolution results
saveRDS(dwls_results, "dwls_deconvolution_results.rds")
write.csv(cell_proportions, "cell_type_proportions.csv")

