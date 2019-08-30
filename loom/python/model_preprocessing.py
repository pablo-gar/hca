'''
Does basic preprocessing of a loom file and returns and H5ad file

Steps of preprocessing:
   - log1p
   - genes with min cells expressed (5% cells)
   - append variable genes logical feature

Notes:
   - Uses the "Gene" in loom features as ids, and makes them unique in case they are not

Usage:
python3 model_preprocessing.py True True True True True in.loom out.h5ad

'''

import scanpy as sc
from anndata import AnnData
from typing import Union, Optional, Tuple, Collection
import sys

MIN_PERCENT_CELLS_GENE_FILTER = 0.05 # Only genes expressed in this fraction of cells
MIN_GENES_FILTER = 10 # Only genes expressed in this fraction of cells
N_VAR_GENES = 2e3 # Only to x variable genes

def main():

    do_log1p = str_to_bool(sys.argv[1])
    do_min_cells_filter = str_to_bool(sys.argv[2])
    do_min_genes_filter = str_to_bool(sys.argv[3])
    do_select_variable_genes = str_to_bool(sys.argv[4])
    do_umap = str_to_bool(sys.argv[5])
    remove_non_variable_genes = str_to_bool(sys.argv[6])

    in_loom_file = sys.argv[-2]
    out_h5ad_file = sys.argv[-1]

    # Read
    print("Reading loom")
    loom = sc.read(in_loom_file)

    # Preprocessing
    loom = preprocessing(loom, do_log1p=do_log1p, do_min_cells_filter=do_min_cells_filter,
                         do_min_genes_filter=do_min_genes_filter,  do_select_variable_genes=do_select_variable_genes,
                         percent_cells=MIN_PERCENT_CELLS_GENE_FILTER, n_var_genes=N_VAR_GENES,
                         min_genes=MIN_GENES_FILTER, remove_non_variable_genes=remove_non_variable_genes)

    if do_umap:
        sc.pp.normalize_total(loom, exclude_highly_expressed=True)
        sc.pp.pca(loom)
        sc.pp.neighbors(loom)
        sc.tl.umap(loom)
        #sc.tl.louvain(loom)

    # Write
    loom.write(out_h5ad_file)

def preprocessing(adata: AnnData,
                  do_log1p: Optional[bool] = True,
                  do_min_cells_filter: Optional[bool] = True,
                  do_min_genes_filter: Optional[bool] = True,
                  do_select_variable_genes: Optional[bool] = True,
                  percent_cells: Optional[float] = 0.05,
                  min_genes: Optional[float] = 10,
                  n_var_genes: Optional[int] = 2e3,
                  remove_non_variable_genes: Optional[bool] = True,
                  flavor: str = 'seurat'
                  ):

    '''
    Wrapper for a series of scanpy preprocessing steps on an annData object

    :param adata: annData
    :param do_log1p: Boolean - If True it will do log1p transformation
    :param do_min_cells_filter: Boolean - If True it will apply gene filter based on no of cells expressing
    :param do_min_genes_filter: Boolean - If True it will apply ce;; filter based on no of genes being expressed
    :param do_select_variable_genes: - If True it will apply select and create variable genes by returning a view of the data
    :param percent_cells: float - frequency of cell expressing genes used for min_cells_filter
    :param min_genes: int - minimum number of expressed genes for a cell to pass the min_genes_filte
    :param n_var_genes: int - select this number of var genes
    :param remove_non_variable_genes: bool -  Only applicable if do_select_variable_genes=True. If True returns AnnData
    only with variable genes, otherwise it returns full object and appends a bool vector to obs metadata
    :param flavor: str -  Flavor for selection of vartiable genes, this is passed to scanpy.pp.highly_variable_genes
    :return: annData
    '''

    adata.var_names_make_unique()

    if do_log1p:
        print("Applying log1p transformation")
        sc.pp.log1p(adata)

    if do_min_cells_filter:
        print("Applying min cells filter: only keep genes if expressed in at least", str(percent_cells * 100), "% cells")
        sc.pp.filter_genes(adata, min_cells=int(percent_cells*adata.n_obs))

    if do_min_genes_filter:
        print("Applying min genes filter: only keep cells if having at least", min_genes, "genes expressed")
        sc.pp.filter_cells(adata, min_genes=min_genes)

    if do_select_variable_genes:
        print("Selecting", str(n_var_genes), "variable genes")
        sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=2000)
        if remove_non_variable_genes:
            highly_variable_genes = adata.var.highly_variable.to_list()
            adata = adata[:, highly_variable_genes]

    return adata


def str_to_bool(x: str):

    if x == "True":
        return True
    elif x == "False":
        return False
    else:
        raise Exception("Only 'True' or 'False' values accepted for options")

if __name__ == "__main__":
    main()
