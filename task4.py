# TASK4.PY     2020/03/02     Elizabeth Chen
# This program performs
#      - scRNA analysis on pancreatic B-cells and cardiomyocytes.
#      - computes average gene expression score for genes related to input keyword.

# ********** 0. IMPORT LIBRARIES ********** #
import scanpy as sc
import numpy as np
import pandas as pd
import csv
# ***************************************** #

sc.settings.verbosity = 3
key_word = 'voltage-gated proton channel'

# ************** 1. LOAD DATA ************* #
# 1a. Extract gene count matrix data for pancreatic B-cells and cardiomyocytes.
adata_pancreas = sc.read_h5ad('tabula-muris-senis-droplet-processed-official-annotations-Pancreas.h5ad')
adata_heartandaorta = sc.read_h5ad('tabula-muris-senis-droplet-processed-official-annotations-Heart_and_Aorta.h5ad')
adata_bcells = adata_pancreas[adata_pancreas.obs['cell_ontology_class']=='pancreatic B cell']
adata_cardiomyocytes = adata_heartandaorta[adata_heartandaorta.obs['cell_ontology_class']=='cardiomyocyte']

# 1b. Combine pancreatic B-cell data with cardiomyocyte data to form one single data object.
adata = adata_bcells.concatenate(adata_cardiomyocytes)
print(adata)
print(adata.var_names)

# 1c. Remove previously created variables
adata.var.drop(labels=['means-0', 'dispersions-0', 'dispersions_norm-0', 'highly_variable-0', 'means-1', 'dispersions-1', 'dispersions_norm-1', 'highly_variable-1'] , axis=1, inplace=True)
print(adata)
# ***************************************** #

# ************ 2. PREPROCESSING *********** #
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata
print(adata)
# ***************************************** #

# *********** 3. MAKE GENE LIST *********** #
data = pd.read_csv('IT_genes.csv')
gene_list = np.array([])

for index, row in data.iterrows():
    if key_word in row['Annotated Term']:
        gene_list = np.append(gene_list, row['Symbol'])
# ***************************************** #

# ************ 4. SCORE GENES ************* #
#print(list(adata.var_names))
print(adata)
sc.tl.score_genes(adata, gene_list)
print(adata)
print(np.shape(adata.obs['score']))
print(adata.obs['score'])

score_b = 0 
score_c = 0
n_b = 0
n_c = 0
for index2, row2 in adata.obs.iterrows():
    if row2['cell_ontology_class'] == 'pancreatic B cell':
        score_b = score_b + row2['score']
        n_b += 1
    else:
        score_c = score_c + row2['score']
        n_c += 1

print('Average score for B cell:', score_b/n_b)
print('Average score for cardiomyocytes:', score_c/n_c)
# ***************************************** #
