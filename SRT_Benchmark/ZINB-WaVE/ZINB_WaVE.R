########################
BiocManager::install("splatter", build_vignettes=TRUE)
BiocManager::install(c("scater"))
library(splatter)
library(zellkonverter)
library(SingleCellExperiment)
library(scater)

# 读取 h5ad 文件
zinb_adata <- readH5AD("D:/2025春夏学期/CMML3/ICA2/scCube/tutorial/demo_data/DLPFC_151507_adata.h5ad")

# 转换为 SingleCellExperiment 对象
zinb_seurat_obj <- as.Seurat(zinb_adata, counts = "X", data = NULL)  # "X" 是 .h5ad 中的主矩阵
zinb_sce <- as.SingleCellExperiment(zinb_seurat_obj)

# 2. 转换 counts 为稠密矩阵
library(Matrix)
counts_matrix <- Matrix::as.matrix(counts(zinb_sce))  # 专用方法，通常保留名称
counts(zinb_sce) <- counts_matrix 

# 3. 拟合参数
params <- splatEstimate(zinb_sce)

# 4. 调整参数（按需修改）
params2 <- setParams(params, 
                    nGenes = 18094,#155
                    batchCells = 4221)#

# 5. 模拟数据
zinb_sim_sce <- splatSimulate(params2, method = "groups")



# 假设原始数据有空间坐标（如 'X' 和 'Y'）
spatial_coords <- colData(zinb_sce)[, c("x", "y")]  # 替换为实际的坐标列名

# 将空间坐标添加到模拟数据
colData(zinb_sim_sce)$x <- spatial_coords$x
colData(zinb_sim_sce)$y <- spatial_coords$y

# 可视化模拟数据的空间分布
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
