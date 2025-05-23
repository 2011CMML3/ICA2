install.packages("devtools")
devtools::install_github("dmcable/spacexr")

library(spacexr)
# load data
reference <- Reference(counts = ref_data$counts, cell_types = ref_data$cell_types)
spatial <- SpatialRNA(counts = spatial_data$counts, coords = spatial_data$coords)

rctd <- create.RCTD(spatial, reference, max_cores = 4)
rctd <- run.RCTD(rctd)

results <- rctd@results
weights <- results$weights
