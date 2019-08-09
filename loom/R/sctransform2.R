# Applies SCtransform to a loom file and creates a new one containing the new normalized layer

library("Seurat")
library("sctransform")
library("loomR")
library("reticulate")

main <- function(cmdArgs = commandArgs(T)) {
    
    loom_in_file <- cmdArgs[1]
    loom_out_file <- cmdArgs[2]
    conda_env <- cmdArgs[3]
    
    use_condaenv(conda_env)
    
    
    #----
    # SCTransform
    
    # Read loom file
    loom_in <- connect(loom_in_file)
    
    # Extract matrix and run sctransform
    mat <- t(loom_in$matrix[,])
    dimnames(mat) <- list(loom_in$row.attrs$Accession[], loom_in$col.attrs$CellID[]) 
    
    # sctransform
    ser <- CreateSeuratObject(mat)
    ser <- SCTransform(ser, return.only.var.genes = F)
    
    loom_in$close_all()
    
    #----
    # Write new loom file with transformation
    
    # Subset from original
    new_mat <- slot(ser$SCT, "counts")
    subset_rows <- which(rownames(mat) %in% rownames(new_mat))
    
    # Creting new loom file with subset dataset
    new_loompy(loom_in_file, loom_out_file, subset_rows - 1)
    
    # Add corrected count layers
    create_add_layer_function()
    py$add_layer(loom_out_file, "SCtransform.normalized", slot(ser$SCT, "counts"))
    py$add_layer(loom_out_file, "SCtransform.normalized.log10", slot(ser$SCT, "data"))
    py$add_layer(loom_out_file, "SCtransform.pearsonResiduals", slot(ser$SCT, "scale.data"))
    
}

new_loompy <- function(original, filename, rows){
    py_run_string(convert = F, code =
        paste0(
'import loompy as lp
from numpy import array
loom_in = lp.connect("', original, '")
loom_out = lp.new("', filename, '")
view = loom_in.view[array([', paste(rows, collapse = ',') ,']),:]
#view = loom_in.view[0:5,:]
loom_out.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)
loom_out.close()
loom_in.close()
'
        )
    )
}

create_add_layer_function <- function() {
    py_run_string(code = "
import loompy as lp
def add_layer(loom_in, layer_name, x):
    loom_in = lp.connect(loom_in)
    loom_in[layer_name] = x
    loom_in.close()
"
                 )
invisible(NULL)
}


main()
