# TASK1.PY     2020/01/16     Elizabeth Chen
# This program performs:
#      - Single-cell RNA sequencing analysis on pancreatic B-cells and cardiomyocytes from
#        the heart.
#      - Removes highly variable genes to successively to reveal more subtle variations in gene profiles of B-cells and cardiomyocytes.

# ********** 0. IMPORT LIBRARIES ********** #
import scanpy as sc
# ***************************************** #

sc.settings.verbosity = 3

# ************* 1. LOAD DATA ************** #
# 1a. Extract gene count matrix data for pancreatic B-cells and cardiomyocytes.
adata_pancreas = sc.read_h5ad('tabula-muris-senis-droplet-processed-official-annotations-Pancreas.h5ad')
adata_heartandaorta = sc.read_h5ad('tabula-muris-senis-droplet-processed-official-annotations-Heart_and_Aorta.h5ad')
adata_bcells = adata_pancreas[adata_pancreas.obs['cell_ontology_class']=='pancreatic B cell']
adata_cardiomyocytes = adata_heartandaorta[adata_heartandaorta.obs['cell_ontology_class']=='cardiomyocyte']

# 1b. Combine pancreatic B-cell data with cardiomyocyte data to form one single data object.
adata = adata_bcells.concatenate(adata_cardiomyocytes)
# ***************************************** #

# *********** 2. PREPROCESSING ************ #
sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
adata.raw = adata
# FIRST HVG FILTERING
adata = sc.pp.filter_genes_dispersion(adata, subset=False, min_disp=0.5, max_disp=None, min_mean=0.0125, max_mean=10, n_bins=20, n_top_genes=None, log=True, copy=True)
sc.pp.log1p(adata)
adata.raw = adata
# ***************************************** #


# **** 7. REMOVE HIGHLY VARIABLE GENES **** #
# FIRST REMOVAL OF HVG
adata = adata[:, adata.var['highly_variable']==False]
print(adata)
# SECOND HVG FILTERING
adata = sc.pp.filter_genes_dispersion(adata, subset=False, min_disp=0.5, max_disp=None, min_mean=0.0125, max_mean=10, n_bins=20, n_top_genes=None, log=True, copy=True)
# SECOND REMOVAL OF HVG
adata = adata[:, adata.var['highly_variable']==False]
print(adata)
# THIRD HVG FILTERING
adata = sc.pp.filter_genes_dispersion(adata, subset=False, min_disp=0.5, max_disp=None, min_mean=0.0125, max_mean=10, n_bins=20, n_top_genes=None, log=True, copy=True)
# THIRD REMOVAL OF HVG
adata = adata[:, adata.var['highly_variable']==False]
print(adata)
# FOURTH HVG FILTERING
#adata = sc.pp.filter_genes_dispersion(adata, subset=False, min_disp=0.5, max_disp=None, min_mean=0.0125, max_mean=10, n_bins=20, n_top_genes=None, log=True, copy=True)
# FOURTH REMOVAL OF HVG
#adata = adata[:, adata.var['highly_variable']==False]
#print(adata)
adata.raw = adata
# ***************************************** #


# ************ 3. PCA ANALYSIS ************ #
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color=['cell_ontology_class'], legend_loc='on data', save='pca.pdf') # Plot
# ***************************************** #

# ********* 4. UMAP VISUALIZATION ********* #
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['cell_ontology_class'], legend_loc='on data', save='umap_by_cell.pdf')
# ***************************************** #

# ******** 5. CLUSTERING ANALYSIS ********* #
# (Optional)??
sc.tl.leiden(adata)
#sc.pl.umap(adata, color=['leiden'], legend_loc='on data', save='umap_by_leiden.pdf')
# ***************************************** #

# ********** 6. FIND MARKER GENES ********* #
sc.tl.rank_genes_groups(adata, 'cell_ontology_class', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# ***************************************** #
