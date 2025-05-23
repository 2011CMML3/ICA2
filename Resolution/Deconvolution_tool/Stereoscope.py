!pip install stereoscope
import stereoscope as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------
# 1. Load Data
# --------------------------

# Generate simulated single-cell data (1000 cells, 200 genes)
sc_data = sc.AnnData(
    X=np.random.negative_binomial(n=5, p=0.1, size=(1000, 200)),  # Simulated UMI counts
    obs={"cell_type": np.random.choice(["Neuron", "Astrocyte", "Microglia"], size=1000)}  # Cell type labels
)

# Generate simulated spatial data (500 spots, 200 genes)
spatial_data = sc.AnnData(
    X=np.random.negative_binomial(n=5, p=0.1, size=(500, 200)),  # Simulated UMI counts
    obsm={"spatial": np.random.rand(500, 2)}  # Spatial coordinates
)

# --------------------------
# 2. Data Preprocessing (Key Steps)
# --------------------------
# Ensure consistent genes between single-cell and spatial data
shared_genes = list(set(sc_data.var_names) & set(spatial_data.var_names))
sc_data = sc_data[:, shared_genes].copy()
spatial_data = spatial_data[:, shared_genes].copy()

# Normalization (Stereoscope requires raw counts, no log transformation)
sc.pp.normalize_total(sc_data, target_sum=1e4)  # Counts per 10,000
sc.pp.normalize_total(spatial_data, target_sum=1e4)

# Convert single-cell data to cell-type × gene pseudo-bulk matrix
cell_type_matrix = pd.crosstab(
    index=sc_data.obs["cell_type"],
    columns=sc_data.var_names,
    values=sc_data.X.toarray().sum(axis=0),  # Assuming X is sparse matrix
    aggfunc="sum"
)
sc_data_pseudo = sc.AnnData(cell_type_matrix)  # Pseudo-bulk data

# --------------------------
# 3. Train Stereoscope Model
# --------------------------
model = st.Stereoscope(
    sc_data=sc_data_pseudo,      # Single-cell pseudo-bulk data (cell types × genes)
    st_data=spatial_data,        # Spatial data
    n_components=10,            # Latent space dimensions (default=10)
    use_gpu=False                # Whether to use GPU acceleration
)
model.fit(
    n_epochs=100,               # Training epochs
    lr=0.01,                    # Learning rate
    batch_size=64               # Batch size
)

# --------------------------
# 4. Predict Spot Cell Type Composition
# --------------------------
results = model.predict()
print(results.head())  # Output cell type proportions per spot (DataFrame)


