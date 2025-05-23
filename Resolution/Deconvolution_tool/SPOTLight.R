devtools::install_github("https://github.com/MarcElosua/SPOTlight")
library(SPOTlight)
sc_data <- readRDS("sc_data.rds")
spatial_data <- readRDS("spatial_data.rds")

# Run SPOTlight
nmf_model <- SPOTlight(
  x = sc_data,
  y = spatial_data,
  groups = sc_data$cell_type
)

# Extract result
cell_proportions <- nmf_model$weights
