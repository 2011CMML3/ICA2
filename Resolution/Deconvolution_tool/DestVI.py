!pip install scvi-tools scanpy

import scvi
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------
# 1. Load Data
# --------------------------
# Load single-cell RNA-seq data (AnnData format)
adata_sc = sc.read("sc_data.h5ad")  

# Load spatial transcriptomics data (AnnData format)
adata_sp = sc.read("spatial_data.h5ad")  

# --------------------------
# 2. Preprocess Data
# --------------------------
# Ensure both datasets use the same genes
common_genes = list(set(adata_sc.var_names) & set(adata_sp.var_names))
adata_sc = adata_sc[:, common_genes].copy()
adata_sp = adata_sp[:, common_genes].copy()

# Normalize single-cell data (counts per 10,000 + log1p)
sc.pp.normalize_total(adata_sc, target_sum=1e4)
sc.pp.log1p(adata_sc)

# Normalize spatial data (same as single-cell)
sc.pp.normalize_total(adata_sp, target_sum=1e4)
sc.pp.log1p(adata_sp)

# --------------------------
# 3. Setup and Train DestVI Model
# --------------------------
# Setup the anndata object for DestVI
scvi.model.DestVI.setup_anndata(
    adata_sc,
    layer="log1p",  
    batch_key=None  
)

# Initialize DestVI model
model = scvi.model.DestVI(
    adata_sc,
    n_latent=10, 
    n_layers=2    
)

# Train the model
model.train(
    max_epochs=100,  
    use_gpu=False    
)

# --------------------------
# 4. Predict Cell Type Proportions
# --------------------------
# Get cell type proportions for spatial data
results = model.predict(
    adata_sp,
    softmax_scale=4.0  # Controls sharpness of cell type assignments
)

# The results object contains:
# - proportions: Cell type proportions per spot
# - concentrations: Cell type-specific expression patterns
cell_proportions = results["proportions"]
print(cell_proportions.head())

# --------------------------
# 5. Visualize Results
# --------------------------
# Add proportions to spatial AnnData object
for ct in cell_proportions.columns:
    adata_sp.obs[f"DestVI_{ct}"] = cell_proportions[ct]

# 6. Save Results
# --------------------------
# Save cell type proportions
cell_proportions.to_csv("destvi_cell_proportions.csv")

# Save model for later use
model.save("destvi_model/")

