!pip install cell2location scanpy numpy pandas matplotlib

import scanpy as sc
import cell2location as c2l
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------
# 1. Load Data
# --------------------------
# Load single-cell RNA-seq data (AnnData format)
# Should contain:
# - Raw counts in .X
# - Cell type annotations in .obs['cell_type']
# - Batch information (if available) in .obs['batch']
adata_sc = sc.read("sc_data.h5ad")

# Load spatial transcriptomics data (AnnData format)
# Should contain:
# - Raw counts in .X
# - Spatial coordinates in .obsm['spatial']
adata_sp = sc.read("spatial_data.h5ad")

# --------------------------
# 2. Preprocess Data
# --------------------------
# Filter genes (keep genes present in both datasets)
shared_genes = list(set(adata_sc.var_names) & set(adata_sp.var_names))
adata_sc = adata_sc[:, shared_genes].copy()
adata_sp = adata_sp[:, shared_genes].copy()

# Normalize single-cell data (CPM + log1p)
sc.pp.normalize_total(adata_sc, target_sum=1e6)
sc.pp.log1p(adata_sc)

# --------------------------
# 3. Model Setup
# --------------------------
# Setup single-cell reference
c2l.models.Cell2location.setup_anndata(
    adata_sc,
    layer=None,  
    batch_key="batch" 
)

# Initialize model
model = c2l.models.Cell2location(
    adata_sc,          
    adata_sp,          
    N_cells_per_location=10,  
    detection_alpha=20        
)

# --------------------------
# 4. Train Model
# --------------------------
model.train(
    max_epochs=100,   
    batch_size=None,  
    train_size=0.9,    
    use_gpu=True      
)

# Plot training loss (check convergence)
plt.plot(model.history["elbo_train"], label="train")
plt.plot(model.history["elbo_test"], label="test")
plt.legend()
plt.ylabel("ELBO Loss")
plt.xlabel("Epochs")
plt.show()

# --------------------------
# 5. Predict Cell Densities
# --------------------------
# Get cell abundance estimates
model.predict(
    use_raw=False,    
    return_samples=False
)

# Export results
results = model.export_results()

# The results dictionary contains:
# - 'q05_cell_abundance_w_sf': 5th percentile of cell abundances
# - 'means_cell_abundance_w_sf': Mean cell abundances
# - 'q95_cell_abundance_w_sf': 95th percentile
cell_densities = results['means_cell_abundance_w_sf']

# --------------------------
# 6. Save Results
# --------------------------
# Save cell abundance estimates
cell_densities.to_csv("cell2location_abundances.csv")

# Save model for later use
model.save("cell2location_model/")

# Save spatial data with results
adata_sp.write("spatial_data_with_cell2location.h5ad")
