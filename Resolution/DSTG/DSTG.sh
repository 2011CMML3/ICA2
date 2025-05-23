################# Convert data
# Python
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('MERFISH_0.06_adata.h5ad')
# extract expression matrix
expression_matrix = adata.X.toarray() 
expression_df = pd.DataFrame(expression_matrix, index=adata.var_names, columns=adata.obs_names)

# extract spatial coordinates
spatial_coords = adata.obs[['x', 'y']]  

# save file as csv file
expression_df.to_csv('MERFISH_0.06_expression_matrix.csv')
spatial_coords.to_csv('MERFISH_0.06_spatial_coordinates.csv')

# R
library(readr)
library(Seurat)

expression_matrix <- read_csv("MERFISH_0.06_expression_matrix.csv", col_names = TRUE)
spatial_coords <- read_csv("MERFISH_0.06_spatial_coordinates.csv", col_names = TRUE)
seurat_obj <- CreateSeuratObject(counts = as.matrix(expression_matrix), meta.data = spatial_coords)

# save as rds file
saveRDS(seurat_obj, "MERFISH_0.06_count.RDS")


################ DSTG
# seqFISH+ MOB internal default gene set
module add python/3.7.9
module add r/4.0.3
module add gcc/4.9.1
module add clang/6.0

cd DSTG
# 转换数据
Rscript convert_data.R ./processed_data/MOB/ref/internal/scRNA_seqfish_count.RDS ./processed_data/MOB/st_mob_seqfish_count.RDS ./processed_data/MOB/ref/internal/scRNA_seqfish_label.RDS
python3 train.py 

# results will be in DSTG_Result/predict_output.csv
# cell type labels can be found in Infor_Data/ST_label/ST_label_1.csv
