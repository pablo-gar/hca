# Creates a loom file from plain texts files: 
#    - A tab-separated or csv file being a gene and cell numeric matrix, and a cell-type file with the same number of lines as number of cells
#
# Assumes that first column is gene names and first row is cell names


import sys, os, loompy
import numpy as np

def main():
    
    # Read args
    mat_file = sys.argv[1]
    cell_features_file = sys.argv[2]
    out_loom_file = sys.argv[3]
    
    delimiter = get_delimiter_from_extension(mat_file)
    
    # Get count matrix, gene names, cell names and 
    count_matrix = read_matrix(mat_file, delimiter)
    cell_names = read_names(mat_file, delimiter, "column")
    gene_names = read_names(mat_file, delimiter, "row")
    
    # Getting features i.e. col and row metadata
    cell_features = read_cell_features(cell_features_file, get_delimiter_from_extension(cell_features_file))
    cell_features["CellID"] = cell_names
    row_features = {"Gene" : gene_names}
    
    # Making sure the number of cells in features matches number of cells in matrix
    ncells = count_matrix.shape[1]
    for key in cell_features:
        if ncells != cell_features[key].shape[0]:
            raise Exception("Number of cells in " + key + "is not the sames as number of cells in matrix")
    
    
    # Make loom file
    loompy.create(out_loom_file, count_matrix, row_attrs = row_features, col_attrs = cell_features,)

def read_matrix(fname, delimiter):
    ncols = get_ncols(fname, delimiter)
    return np.loadtxt(fname = fname, delimiter = delimiter, skiprows = 1, usecols = range(1, ncols))

def read_names(fname, delimiter, dimension):
    '''
    Reads row names (dimension = "row") or column names (dimension = "column")
    of a file into a numpy array 
    '''
    
    with open(fname, "r") as f:
        if dimension == "column":
            return np.array(f.readline().rstrip().split(delimiter)[1:])
        elif dimension == "row":
            return np.loadtxt(fname = fname, dtype = str, delimiter = delimiter, skiprows = 1, usecols = 0)
        else:
            raise('dimension has to be "row" or "column"')
 
def get_ncols(fname, delimiter):
    with open(fname, "r") as f:
        return len(f.readline().split(delimiter))
    
def read_cell_features(fname, delimiter):
    '''
    Reads from text file a table with each row corresponding to individual cell features
    Number of rows should be number of features
    Number of columns should be number of cells + 1
    First  column should be name of features
    
    Right now only string features are supported
    
    returns a dictionary with keys being cell feature names and values numpy arrays
    '''
    
    x = dict() 
    with open(fname, "r") as f:
        for line in f:
            line = line.rstrip().split(delimiter)
            x[line[0]] = np.array(line[1:])
    
    return x

def get_delimiter_from_extension(fname):
    ext = os.path.splitext(fname)[1]
    if ext == ".tsv":
        delimiter = "\t"
    elif ext == ".csv":
        delimiter = ","
    else:
        raise("Only .tsv or .csv matrices supported")
    
    return delimiter
            
if __name__ == "__main__":
    main()
    
   
