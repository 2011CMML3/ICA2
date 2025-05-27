setwd('D:/CMML3/ICA2/R')
getwd()
library(Seurat)
library(SingleCellExperiment)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("zellkonverter")

# Read the file
library(zellkonverter)
# BiocManager::install("mkl-2024.2.2")
sce <- readH5AD("D:/CMML3/ICA2/scCube/tutorial/demo_data/DLPFC_151507_adata.h5ad") # MERFISH_0.06_adata.h5ad
class(sce)
print(sce)
# Convert to seurat object
library(Seurat)
seurat_obj <- as.Seurat(sce, counts = "X", data = NULL)  
sce <- as.SingleCellExperiment(seurat_obj)
sce2 <- sce



library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(viridis)
theme_set(theme_classic())
genes <- c('GFAP','HPCAL1','HOPX','NEFH','PCP4','KRT17','MOBP') #'Gad1','Mbp','Nnat','Aqp4','Slc17a6','Fn1','Pdgfra','Selplg','Myh11':for MERFISH data
example_sce <- sce2[which(rownames(sce)%in%genes), ]
# set.seed(123)
# simulate the data
example_simu <- scdesign3(
  sce = example_sce,
  assay_use = "counts",
  celltype = "Cell_type",
  pseudotime = NULL,
  spatial = c("x", "y"),#spatial1,2-> x,y
  other_covariates = NULL,
  mu_formula = "s(x, y, bs = 'gp', k= 400)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "pbmcapply"
)

# convert to singlecellexperiment
simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count), colData = example_simu$new_covariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))

# For visualization, extract gene expression profiles and spatial patterns
VISIUM_dat_test <- data.frame(t(log1p(counts(example_sce)))) %>% as_tibble() %>% dplyr::mutate(X = colData(example_sce)$x, Y = colData(example_sce)$y) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "Reference")
VISIUM_dat_scDesign3 <- data.frame(t(log1p(counts(simu_sce)))) %>% as_tibble() %>% dplyr::mutate(X = colData(simu_sce)$x, Y = colData(simu_sce)$y) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "scDesign3")
VISIUM_dat <- bind_rows(VISIUM_dat_test, VISIUM_dat_scDesign3) %>% dplyr::mutate(Method = factor(Method, levels = c("Reference", "scDesign3")))

VISIUM_dat %>% filter(Gene %in% genes) %>% ggplot(aes(x = X, y = Y, color = Expression)) + geom_point(size = 0.5) + scale_colour_gradientn(colors = viridis_pal(option = 'viridis')(10), limits=c(0, 4)) +  theme_classic()  + facet_grid(Method ~ Gene)#+ facet_wrap(~ Gene, ncol = 4)
