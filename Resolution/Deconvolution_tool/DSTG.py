!pip install DSTG

from DSTG import DSTG
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------
# 1. Load Data
# --------------------------
# Load single-cell RNA-seq data (AnnData format)
# Should contain:
# - Normalized counts in .X
# - Cell type annotations in .obs['cell_type']
sc_data = sc.read("sc_data.h5ad")

# Load spatial transcriptomics data (AnnData format)
# Should contain:
# - Normalized counts in .X
# - Spatial coordinates in .obsm['spatial']
spatial_data = sc.read("spatial_data.h5ad")

# --------------------------
# 2. Preprocess Data
# --------------------------
# Ensure both datasets use the same genes
shared_genes = list(set(sc_data.var_names) & set(spatial_data.var_names))
sc_data = sc_data[:, shared_genes].copy()
spatial_data = spatial_data[:, shared_genes].copy()

# Normalize data (CPM + log1p recommended)
sc.pp.normalize_total(sc_data, target_sum=1e4)
sc.pp.log1p(sc_data)
sc.pp.normalize_total(spatial_data, target_sum=1e4)
sc.pp.log1p(spatial_data)

# --------------------------
# 3. Train DSTG Model
# --------------------------
model = DSTG(
    sc_data=sc_data,              # Single-cell reference
    sc_label='cell_type',         # Column name for cell type labels
    spatial_data=spatial_data,    # Spatial transcriptomics data
    use_graph=True,              # Whether to use graph structure (recommended)
    n_neighbors=15,              # Number of spatial neighbors for graph
    use_gpu=False                 # GPU acceleration if available
)

model.train(
    n_epochs=300,                # Training iterations
    lr=0.001,                   # Learning rate
    batch_size=64               # Batch size
)

# --------------------------
# 4. Predict Cell Type Proportions
# --------------------------
results = model.predict()
print(results.head())  # DataFrame with spot × cell_type proportions


# --------------------------
# 5. Save Results
# --------------------------
# Save cell type proportions
results.to_csv("dstg_cell_proportions.csv")

# Save model parameters
model.save("dstg_model/")

