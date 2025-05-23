# Install required packages
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("https://github.com/MarcElosua/SPOTlight")

# Load libraries
library(Seurat)
library(SPOTlight)
library(ggplot2)

# --------------------------
# 1. Load Data
# --------------------------
# Load single-cell RNA-seq data (10X format)
sc_data <- Read10X(data.dir = "sc_data/")
sc_meta <- read.csv("sc_data/metadata.csv")  # Should contain 'cell_type' column
sc_seurat <- CreateSeuratObject(counts = sc_data, meta.data = sc_meta)

# Load spatial transcriptomics data (10X Visium format)
spatial_data <- Read10X(data.dir = "spatial_data/")
spatial_seurat <- CreateSeuratObject(counts = spatial_data)

# Add spatial coordinates (example for Visium)
coordinates <- read.csv("spatial_data/tissue_positions.csv")
rownames(coordinates) <- coordinates$barcode
spatial_seurat@images$slice1 <- new(
  Class = "VisiumV1",
  coordinates = coordinates[,c("x", "y")]
)

# --------------------------
# 2. Preprocess Data
# --------------------------
# Single-cell preprocessing
sc_seurat <- NormalizeData(sc_seurat)
sc_seurat <- FindVariableFeatures(sc_seurat)
sc_seurat <- ScaleData(sc_seurat)

# Spatial data preprocessing
spatial_seurat <- NormalizeData(spatial_seurat)
spatial_seurat <- ScaleData(spatial_seurat)

# Ensure common features
features <- intersect(VariableFeatures(sc_seurat), rownames(spatial_seurat))
sc_seurat <- subset(sc_seurat, features = features)
spatial_seurat <- subset(spatial_seurat, features = features)

# --------------------------
# 3. Run SPOTlight Deconvolution
# --------------------------
deconv_results <- SPOTlight(
  x = GetAssayData(sc_seurat, slot = "data"),  # Normalized scRNA-seq data
  y = GetAssayData(spatial_seurat, slot = "data"),  # Normalized spatial data
  clusters = sc_seurat$cell_type,  # Cell type labels from metadata
  hvg = features,  # Use common highly variable genes
  min_cont = 0.01  # Minimum contribution threshold
)

# --------------------------
# 4. Extract and Visualize Results
# --------------------------
# Get cell type proportions matrix
cell_fractions <- deconv_results$mat
head(cell_fractions)

# Add results to spatial Seurat object
spatial_seurat[["SPOTlight"]] <- CreateAssayObject(data = t(cell_fractions))

# Visualize T cell proportions
SpatialFeaturePlot(
  spatial_seurat,
  features = "T_cell",  # Replace with your cell type name
  pt.size.factor = 1.5,
  alpha = c(0.1, 1)
) +
  scale_fill_gradientn(
    colours = c("lightgrey", "navyblue"),
    limits = c(0, 1),
    name = "Proportion"
  ) +
  ggtitle("T Cell Spatial Distribution")

# --------------------------
# 5. Save Results
# --------------------------
# Save deconvolution results
saveRDS(deconv_results, "spotlight_deconvolution_results.rds")
write.csv(cell_fractions, "cell_type_proportions.csv")

# Save spatial data with deconvolution results
saveRDS(spatial_seurat, "spatial_data_with_deconvolution.rds")
