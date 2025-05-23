!pip install tangram-sc

import scanpy as sc
import tangram as tg

# load data
adata_sc = sc.read("sc_data.h5ad")  # single cell data
adata_sp = sc.read("spatial_data.h5ad")  # spatial data

# map cell types
tg.pp_adatas(adata_sc, adata_sp, genes=None)
ad_map = tg.map_cells_to_space(adata_sc, adata_sp)

# deconvolution results
tg.project_cell_annotations(ad_map, adata_sp, annotation="cell_type")
adata_sp.obs.head()  # the proportions of cell types
