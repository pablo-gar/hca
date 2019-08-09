# Applies SCtransform to a loom file and creates a new one containing the new normalized layer

library("Seurat")
library("sctransform")
library("loomR")

main <- function(cmdArgs = commandArgs(T)) {
    
    loom_in_file <- cmdArgs[1]
    loom_out_file <- cmdArgs[2]
    
    # Read loom file
    loom_in <- connect(loom_in_file)
    
    #----
    # SCTransform
    
    # Extract matrix and run sctransform
    mat <- t(loom_in$matrix[,])
    dimnames(mat) <- list(loom_in$row.attrs$Accession[], loom_in$col.attrs$CellID[]) 
    
    # sctransform
    ser <- CreateSeuratObject(mat)
    ser <- SCTransform(ser, return.only.var.genes = F)
    
    #----
    # Write new loom file with transformation
    
    # Subset from original
    subset_rows <- which(rownames(mat) %in% rownames(slot(ser$SCT, "counts")))
    subset(loom_in, n = subset_rows, filename = loom_out_file)
    loom_out <- connect(loom_out_file, "r+")
    
    # Add corrected count layers
    
    row_ids <- loom_out$row.attrs$Accession[]
    col_ids <- loom_out$col.attrs$CellID[]
    
    layers <- list(SCTransform.counts = slot(ser$SCT, "counts")[row_ids, col_ids],
                   SCTransform.log10Counts = slot(ser$SCT, "data")[row_ids, col_ids],
                   SCTransform.pearsonResiduals = slot(ser$SCT, "scale.data")[row_ids, col_ids]
    )
    
    loom_out <- loom_out$add.layer(layers)
    
    loom_out$close_all()
    loom_in$close_all()
}


main()