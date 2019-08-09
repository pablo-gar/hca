from scvi.models import SCANVI
from scvi.dataset.anndataset import AnnDatasetFromAnnData
from scvi.inference import UnsupervisedTrainer, JointSemiSupervisedTrainer, SemiSupervisedTrainer
import scanpy as sc
import numpy as np
import torch

in_loom_train = "/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/paper.loom"
#in_loom_train = "/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/paper_.loom"
in_loom_test = "/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/hca_cellTypes.loom"

#---------------------------
# Example of pre-processing

adata = sc.read_loom(in_loom_train)
adata_test = sc.read_loom(in_loom_test)

adata.var_names_make_unique()
adata_test.var_names_make_unique()


# PRE-PROCESS
# First find do log1p
sc.pp.filter_genes(adata, min_cells=int(0.05 * adata.shape[0]))
sc.pp.log1p(adata)
sc.pp.log1p(adata_test)
# Then find variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
# Label test data
adata.obs["scanvi_test"] = False
adata_test.obs["scanvi_test"] = True

# SELECT SAME GENES
genes = adata[:,adata.var.highly_variable.tolist()].var_names
genes_shared = [i in adata_test.var_names.to_list() for i in genes]
genes = genes[genes_shared]

adata = adata[:, genes.to_list()]
adata_test = adata_test[:, genes.to_list()]

adata_merged = adata.concatenate(adata_test)

# SCANVI
adata_scanvi = AnnDatasetFromAnnData(adata_merged)
scanvi = SCANVI(adata_scanvi.nb_genes, adata_scanvi.n_batches, adata_scanvi.n_labels)
trainer_scanvi = SemiSupervisedTrainer(scanvi, adata_scanvi, frequency=5)

n_epochs = 200
trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=adata_scanvi.)