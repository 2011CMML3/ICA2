!pip install tangram-sc scanpy numpy pandas matplotlib
import scanpy as sc
import tangram as tg
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------
# 1. Load Data
# --------------------------
# Load single-cell RNA-seq data (AnnData format)
# Should contain:
# - Normalized expression in .X
# - Cell type annotations in .obs['cell_type']
adata_sc = sc.read("sc_data.h5ad")

# Load spatial transcriptomics data (AnnData format)
# Should contain:
# - Normalized expression in .X
# - Spatial coordinates in .obsm['spatial']
adata_sp = sc.read("spatial_data.h5ad")

# --------------------------
# 2. Preprocess Data
# --------------------------
# Ensure both datasets use the same genes
tg.pp_adatas(
    adata_sc, 
    adata_sp,
    genes=None, 
)

# --------------------------
# 3. Map Cells to Space
# --------------------------
# Train Tangram alignment model
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_sp,
    mode='clusters',  # Alternative: 'cells' for single-cell resolution
    cluster_label='cell_type',  # Column name with cell annotations
    density_prior='rna_count_based',  # Weight by RNA abundance
    num_epochs=1000,  # Training iterations
    device='cuda' if tg.utils.is_cuda_available() else 'cpu'
)

# --------------------------
# 4. Project Cell Type Annotations
# --------------------------
# Transfer cell type labels to spatial data
tg.project_cell_annotations(
    ad_map, 
    adata_sp,
    annotation='cell_type'  # Column name in adata_sc.obs
)

# View results (cell type proportions per spot)
print(adata_sp.obs.head())


# 5. Save Results
# --------------------------
# Save spatial data with deconvolution results
adata_sp.write('spatial_data_with_tangram.h5ad')

# 
spatial_coords = pd.DataFrame(adata_sp.obsm['spatial'], 
                             columns=['x', 'y'], 
                             index=adata_sp.obs_names)
result_df = pd.concat([spatial_coords, adata_sp.obs], axis=1)
result_df.to_csv('Tangram_celltype_results.csv')
