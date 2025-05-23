devtools::install_github("RubD/Giotto")
library(Giotto)

# load data
sc_data <- readRDS("sc_data.rds")  # single cell data
spatial_data <- readRDS("spatial_data.rds")  # spatial data

# Run spatialDWLS
dwls_results <- runDWLSDeconv(
  spatial_data,
  sc_data,
  cell_type_column = "cell_type"
)

# Extract result
cell_proportions <- dwls_results$weights
