library(Seurat)
library(SeuratObject)
library(dplyr)


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
data_dir <- args[1]

# create a Seurat object
seu_object <- Seurat::Load10X_Spatial(data.dir = data_dir)

# write to a seperate
seu_meta <- seu_object@meta.data

# read in GSVA table
GSVA_paths <- readRDS("per_spot_scores.RDS")
GSVA_paths$Barcode = gsub(".bam","",GSVA_paths$Barcode)
geneset_ids = colnames(GSVA_paths)[-1]


# copy the rownames column (the 10x barcodes) to a different column in the meta_data. This is so we can check that the join was done correctly.
seu_meta_new <- seu_meta %>%  mutate(Barcode = rownames(.))

## now join the pathways df to the extracted meta-data df
meta_data_joined <- left_join(x = seu_meta_new, y = GSVA_paths, by = "Barcode", keep=T)

for (geneset in geneset_ids){
  seu_object <- AddMetaData(object = seu_object, metadata = meta_data_joined, col.name = geneset)
  
  
}
     
saveRDS(seu_object,"Seurat_output.RDS")
SpatialPlot(object = seu_object, features = geneset_ids, pt.size.factor = 2,combine=FALSE)
dev.off() #SpatialPlot will write to a PDF, but it won't save and close the file until you turn off the plotting device.
     
     
