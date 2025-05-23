install.packages("Seurat")
library(Seurat)
library(SPOTlight)

# load data
sc_data <- Read10X("sc_data/")
spatial_data <- Read10X("spatial_data/")

# Run SPOTlight
deconv_results <- SPOTlight(
  sc_data,
  spatial_data,
  clusters = "cell_type"
)

# Extract result
cell_fractions <- deconv_results$mat
