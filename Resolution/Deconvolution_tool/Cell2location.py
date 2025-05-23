!pip install cell2location

import scanpy as sc
import cell2location as c2l

# load data
adata_sc = sc.read("sc_data.h5ad")
adata_sp = sc.read("spatial_data.h5ad")

# train the models using default parameters
c2l.models.Cell2location.setup_anndata(adata_sc)
model = c2l.models.Cell2location(adata_sc, adata_sp)
model.train(max_epochs=100)

# predict cell density
model.predict()
results = model.export_results()
