!pip install scvi-tools

import scvi
import scanpy as sc

# load data
adata_sc = sc.read("sc_data.h5ad")
adata_sp = sc.read("spatial_data.h5ad")

# Train DestVI model
scvi.model.DestVI.setup_anndata(adata_sc)
model = scvi.model.DestVI(adata_sc)
model.train(max_epochs=100)

# predict spot composition
results = model.predict(adata_sp)
