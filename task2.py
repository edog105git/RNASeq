# TASK2.PY     2020.02.02     Elizabeth Chen
# This program performs:
#      - scRNA analysis on pancreatic B-cells and cardiomyocytes.
#      - Keeps data only from genes with select prefixes (e.g. CACNA) and performs analysis (see Section 3 of code).

# ********** 0. IMPORT LIBRARIES ********** #
import scanpy as sc
import numpy as np
# ***************************************** #

sc.settings.verbosity = 3

# ************** 1. LOAD DATA ************* #
# 1a. Extract gene count matrix data for pancreatic B-cells and cardiomyocytes.
adata_pancreas = sc.read_h5ad('tabula-muris-senis-droplet-processed-official-annotations-Pancreas.h5ad')
adata_heartandaorta = sc.read_h5ad('tabula-muris-senis-droplet-processed-official-annotations-Heart_and_Aorta.h5ad')
adata_bcells = adata_pancreas[adata_pancreas.obs['cell_ontology_class']=='pancreatic B cell']
adata_cardiomyocytes = adata_heartandaorta[adata_heartandaorta.obs['cell_ontology_class']=='cardiomyocyte']

# 1b. Combine pancreatic B-cell data with cardiomyocyte data to form one single data object.
adata = adata_bcells.concatenate(adata_cardiomyocytes)
print(adata)
print(adata.var)

# 1c. Remove previously created variables
adata.var.drop(labels=['n_cells', 'means-0', 'dispersions-0', 'dispersions_norm-0', 'highly_variable-0', 'means-1', 'dispersions-1', 'dispersions_norm-1', 'highly_variable-1'] , axis=1, inplace=True)
print(adata)
# ***************************************** #



# ************ 2. PREPROCESSING *********** #
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata
print(adata)
# ***************************************** #

# ** 3. SELECT FOR GENES WITH SPECIFIC PREFIX ** #
ca_channel = []
genes_present = list(adata.var_names)

for id in genes_present:
    ca_channel.append(id.lower().startswith('cacn'))

adata.var['ca_channel'] = ca_channel
print(adata.var)
adata = adata[:, adata.var['ca_channel']==True]
adata.raw = adata
print(adata)
# ********************************************** #        


# ********** 4. PCA ANALYSIS ********** #
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color=['cell_ontology_class'], legend_loc='on data')
# ************************************* #

# ******* 5. UMAP VISUALIZATION ******* #
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['cell_ontology_class'], legend_loc='on data')
# ************************************* #

# ********** 6. FIND MARKERS ********** #
sc.tl.rank_genes_groups(adata, 'cell_ontology_class', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False)
# ************************************* #
