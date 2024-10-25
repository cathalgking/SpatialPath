library(Seurat)
library(SeuratObject)
library(dplyr)

# the directory containing the SpaceRanger "/outs/" folder
data_dir <- "./outs/"

# create a Seurat object
seu_object <- Seurat::Load10X_Spatial(data.dir = data_dir)

# write to a seperate
seu_meta <- seu_object@meta.data

# read in GSVA table
GSVA_paths <- readRDS("./outs/per_spot_scores.RDS")
GSVA_paths$Barcode = gsub(".bam","",GSVA_paths$Barcode)
#GSVA_paths <- GSVA_paths[,1:2] # only keep 1 pathway for simplicity


# copy the rownames column (the 10x barcodes) to a different column in the meta_data. This is so we can check that the join was done correctly.
seu_meta_new <- seu_meta %>%  mutate(Barcode = rownames(.))

## now join the pathways df to the extracted meta-data df
meta_data_joined <- left_join(x = seu_meta_new, y = GSVA_paths, by = "Barcode", keep=T)

## add the meta-data column with the matched barcodes to the seurat object. each pathway 1-by-1
seu_object <- AddMetaData(object = seu_object, metadata = meta_data_joined, col.name = "GOBP_INNATE_IMMUNE_RESPONSE")
seu_object <- AddMetaData(object = seu_object, metadata = meta_data_joined, col.name = "GOBP_ADAPTIVE_IMMUNE_RESPONSE")
seu_object <- AddMetaData(object = seu_object, metadata = meta_data_joined, col.name = "GOBP_HORMONE_METABOLIC_PROCESS")
seu_object <- AddMetaData(object = seu_object, metadata = meta_data_joined, col.name = "GOBP_MITOTIC_CELL_CYCLE")
seu_object <- AddMetaData(object = seu_object, metadata = meta_data_joined, col.name = "GOBP_ATP_METABOLIC_PROCESS")
seu_object <- AddMetaData(object = seu_object, metadata = meta_data_joined, col.name = "GOBP_PROTEIN_FOLDING")


saveRDS(seu_object,"Seurat_output.RDS")
SpatialPlot(object = seu_object, features = c("GOBP_INNATE_IMMUNE_RESPONSE",
                                              "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
                                              "GOBP_HORMONE_METABOLIC_PROCESS",
                                              "GOBP_ATP_METABOLIC_PROCESS",
                                              "GOBP_PROTEIN_FOLDING",
                                              "GOBP_MITOTIC_CELL_CYCLE"
                                              ), pt.size.factor = 2,combine=FALSE)
dev.off() #SpatialPlot will write to a PDF, but it won't save and close the file until you turn off the plotting device.
