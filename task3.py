# TASK3.PY     2020/02/06     Elizabeth Chen
# This program performs
#      - scRNA analysis on pancreatic B-cells and cardiomyocytes.
#      - Only performs analysis on genes from an input gene list csv file.
#      - Provide name of input gene list csv file in Section 3.

# ********** 0. IMPORT LIBRARIES ********** #
import scanpy as sc
import numpy as np
import pandas as pd
import csv
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
data = pd.read_csv('IT_genes_v2.csv')
gene_list = list(data['Symbol'])

with open('gene_list.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    for i in range(len(gene_list)):
        writer.writerow([gene_list[i]])
# ***************************************** #

# ** 4. SELECT FOR GENES FROM INPUT LIST ** #
genes_present = list(adata.var_names)
channel_genes = []

for id in genes_present:
    present = False
    if id in gene_list:
#        print(id)
        present = True
    channel_genes.append(present)
    
adata.var['ion channel'] = channel_genes
print(adata.var)
adata = adata[:, adata.var['ion channel']==True]
adata.raw = adata
print(adata)
# ***************************************** #

# ************ 5. PCA ANALYSIS ************ #
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color=['cell_ontology_class'], legend_loc='on data')
# ***************************************** #

# ********* 6. UMAP VISUALIZATION ********* #
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['cell_ontology_class'], legend_loc='on data')
# ***************************************** #

# ************ 7. RANK MARKERS ************ #
sc.tl.rank_genes_groups(adata, 'cell_ontology_class', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False)
# ***************************************** #
