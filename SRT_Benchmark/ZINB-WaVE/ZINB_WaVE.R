########################
BiocManager::install("splatter", build_vignettes=TRUE)
BiocManager::install(c("scater"))
library(splatter)
library(zellkonverter)
library(SingleCellExperiment)
library(scater)

# Read the file
zinb_adata <- readH5AD("D:/2025春夏学期/CMML3/ICA2/scCube/tutorial/demo_data/DLPFC_151507_adata.h5ad")

# Convert SingleCellExperiment object
zinb_seurat_obj <- as.Seurat(zinb_adata, counts = "X", data = NULL)  # "X" 是 .h5ad 中的主矩阵
zinb_sce <- as.SingleCellExperiment(zinb_seurat_obj)

# 
library(Matrix)
counts_matrix <- Matrix::as.matrix(counts(zinb_sce))  # 专用方法，通常保留名称
counts(zinb_sce) <- counts_matrix 

# simulate parameters
params <- splatEstimate(zinb_sce)

# adjust the parameters
params2 <- setParams(params, 
                    nGenes = 18094,#155
                    batchCells = 4221)#5232

# sumulate the data
zinb_sim_sce <- splatSimulate(params2, method = "groups")



# adjust the coordinates
spatial_coords <- colData(zinb_sce)[, c("x", "y")]  

# add coordinates to the simulated data
colData(zinb_sim_sce)$x <- spatial_coords$x
colData(zinb_sim_sce)$y <- spatial_coords$y

# visualize the distribution of the spatial data
library(ggplot2)
ggplot(as.data.frame(colData(zinb_sim_sce)), aes(x, y)) +
  geom_point(size=1, alpha=0.6) +
  theme_bw()

#
VISIUM_dat_test <- data.frame(t(log1p(counts(zinb_sce)))) %>% as_tibble() %>% dplyr::mutate(X = colData(zinb_sce)$x, Y = colData(zinb_sce)$y) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "Reference")

VISIUM_dat_zinb <- data.frame(t(log1p(counts(zinb_sim_sce)))) %>% as_tibble() %>% dplyr::mutate(X = colData(zinb_sim_sce)$x, Y = colData(zinb_sim_sce)$y) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "ZINB-WaVE")
max(VISIUM_dat_zinb$Expression)#4.09#3.9
VISIUM_dat_zinb$Gene <- VISIUM_dat_test$Gene 
VISIUM_dat <- bind_rows(VISIUM_dat_test, VISIUM_dat_zinb) %>% dplyr::mutate(Method = factor(Method, levels = c("Reference",  "ZINB-WaVE")))
VISIUM_dat %>% filter(Gene %in% genes) %>% ggplot(aes(x = X, y = Y, color = Expression)) + geom_point(size = 0.5) + scale_colour_gradientn(colors = viridis_pal(option = 'D')(10), limits=c(0, 4)) + theme_classic()+ facet_grid(Method ~ Gene )
