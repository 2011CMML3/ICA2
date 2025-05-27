if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("dmcable/spacexr")
library(spacexr)
library(Matrix) 
library(ggplot2)

# --------------------------
# 1. Prepare Reference Data
# --------------------------
# Load single-cell reference data
# Should contain:
# - counts: gene x cell sparse matrix (dgCMatrix)
# - cell_types: named vector of cell type assignments
ref_data <- readRDS("sc_reference.rds")

# Create Reference object
reference <- Reference(
  counts = ref_data$counts,        
  cell_types = ref_data$cell_types,
  nUMI = colSums(ref_data$counts)  
)

# --------------------------
# 2. Prepare Spatial Data
# --------------------------
# Load spatial transcriptomics data
# Should contain:
# - counts: gene x spot sparse matrix
# - coords: data.frame with x,y coordinates
spatial_data <- readRDS("spatial_data.rds")

# Create SpatialRNA object
spatial <- SpatialRNA(
  counts = spatial_data$counts,  
  coords = spatial_data$coords, 
  nUMI = colSums(spatial_data$counts)  
)

# --------------------------
# 3. Run RCTD Deconvolution
# --------------------------
# Create RCTD object
rctd <- create.RCTD(
  spatial = spatial,
  reference = reference,
  max_cores = 4,          # Parallel processing
  CELL_MIN_INSTANCE = 10,  # Minimum cells per type
  gene_cutoff = 0.000125,  # Gene expression threshold
  fc_cutoff = 0.5,         # Fold-change threshold
  UMI_min = 100            # Minimum UMIs per spot
)

# Run deconvolution
rctd <- run.RCTD(rctd)

# --------------------------
# 4. Extract and Analyze Results
# --------------------------
# Get full results object
results <- rctd@results

# Extract cell type proportions (weights matrix)
weights <- results$weights
head(weights)

# Get cell type assignments per spot
cell_type_df <- results$results_df
head(cell_type_df)

# 5. Save Results
# --------------------------
# Save weights matrix
write.csv(weights, "rctd_cell_type_weights.csv")



# Save spatial data with deconvolution results
spatial_data$rctd_weights <- weights
saveRDS(spatial_data, "spatial_data_with_rctd.rds")
