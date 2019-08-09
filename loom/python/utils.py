# Custome utils to work with loom objects

import loompy as lp
import scanpy as sp
from itertools import compress

def get_matrix(x, rows_bool = None, cols_bool = None):
    
    '''
    Returns the count matrix using a selection object to slice a numpy array
    
    x - loom object
    selection - obj to slice the 2d numpy array
    
    '''
    
    if (rows_bool is not None) and (cols_bool is not None):
        mat = x[get_true_indices(rows_bool), get_true_indices(cols_bool)]
    elif (rows_bool is not None) and (cols_bool is None):
        mat = x[get_true_indices(rows_bool),:]
    elif (rows_bool is None) and (cols_bool is not None):
        mat = x[:,get_true_indices(cols_bool)]
    elif (rows_bool is None) and (cols_bool is None):
        mat = x[:,:]
    
    return mat


def get_true_indices(t):
    '''
    Gets the indices whhere is T of a boolean list
    '''
    
    return list(compress(range(len(t)), t))


def find_variable_genes(filename, var_names = 'Accession'):
    
    '''
    Appends a boolean row attribute to a loom file indicating if the gene is
    amongst the most variable ones
    '''
    
    x = sp.read_loom(filename, var_names=var_names)

    sp.pp.log1p(x)
    sp.pp.highly_variable_genes(x)
