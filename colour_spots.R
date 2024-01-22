## colour spots based on a selection of barcodes

library(Seurat)
library(SpatialExperiment)
library(dplyr)
library(tibble)
#library(ggspavis)
#library(Voyager)


#####
## Seurat object
A1_seu <- readRDS("/Users/cathal.king/Documents/Projects/DL/Visium_Mouse/Seurat/Visium_mouse_Re_seurat_A1.rds")
SpatialPlot(object = A1_seu)


# path must contain a 10x "/outs/" folder
# sample <- file.path("/Users/cathal.king/Documents/Projects/DL/Visium_Mouse/SpaceRanger_trial_grey/outs")
# # read in
# A1_spe <- read10xVisium(samples = sample, type = "sparse", images = "lowres")
# A1_spe
# # View image in SPE object
# plotVisium(A1_spe)


###########################################################
## pathway df
A_pathways <- read.csv(file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/GSVA_results/A1_Per_spot_activity_scores.csv", row.names = 1)
# colnames of the pathways df are the barcodes
# rownames of the pathways df are the pathways

# test
# Remove the "_A1" from the end of each column name
colnames(A_pathways) <- sub("_A1$", "", colnames(A_pathways))

A_pathways_t <- as.data.frame(t(A_pathways))

# copy the rownames column (the 10x barcodes) to a different column in the pathways df
A_pathways_t <- A_pathways_t %>%
  mutate(Barcode = rownames(.))

###########################################################

###########################################################
### if using a select list of barcodes from 10x Loupe browser
#barcodes_to_add <- read.table(file = "/Users/cathal.king/Documents/Projects/DL/spatial_plots/barcodes2.txt", header = F, col.names = "barcodes")
barcodes_to_add <- read.csv(file = "/Users/cathal.king/Documents/Projects/DL/spatial_plots/A1_tumour_barcodes.csv", header = T, sep = ",")
# drop the other column
barcodes_to_add <- barcodes_to_add %>% select(-A1_tumour)

# test
# Remove the "-1" from the end of each column name
colnames(A1_seu) <- sub("-1$", "", colnames(A1_seu))
###########################################################


# extract seu metadata
meta_data_A1 <- A1_seu@meta.data

# copy the rownames column (the 10x barcodes) to a different column in the meta_data
meta_data_A1_new <- meta_data_A1 %>%
  mutate(Barcode = rownames(.))

#############
## now join the pathways df to the extracted meta-data df
meta_data_A1_joined <- left_join(x = meta_data_A1_new, y = A_pathways_t2, by = "Barcode", keep=T)


### add that new meta data back into the Seurat object

## add the meta-data column with the matched barcodes to the seurat object. each pathway 1-by-1
A1_seu <- AddMetaData(object = A1_seu, metadata = meta_data_A1_joined, col.name = "GOBP_LYMPHOCYTE_ACTIVATION")

A1_seu <- AddMetaData(object = A1_seu, metadata = meta_data_A1_joined, col.name = "GOBP_B_CELL_ACTIVATION")
A1_seu <- AddMetaData(object = A1_seu, metadata = meta_data_A1_joined, col.name = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION")
A1_seu <- AddMetaData(object = A1_seu, metadata = meta_data_A1_joined, col.name = "GOBP_HORMONE_METABOLIC_PROCESS")




















################################# the below if for if using the 10x Loupe browser file

# copy the rownames column (the 10x barcodes) to a different column in the meta_data
meta_data_A1_new <- meta_data_A1 %>%
  mutate(Barcode = rownames(.))

##### add barcode column to that extracted meta-data
meta_data_A1_joined <- left_join(x = meta_data_A1_new, y = barcodes_to_add, by = "Barcode", keep=T)

####################### convert NA's to "other"
############ or this ## this will definitely convert any NA entries to "other"
data_modified_no_NA <- meta_data_A1_joined %>%
  mutate(category = ifelse(is.na(Barcode.y), "other", Barcode.y))

############ this definitely converts all entries that are not "other" to "tumour"
data_modified_no_NA_incl_Tumour <- data_modified_no_NA %>%
  mutate(category_new = ifelse(category != "other", "tumour", category))


##### combine it into 1 command?
combined_data <- meta_data_A1_joined %>%
  mutate(category_to_colour_by = ifelse(is.na(Barcode.y), "other", 
                               ifelse(Barcode.y != "other", "tumour", Barcode.y)))

## add the meta-data column with the matched barcodes to the seurat object
A1_seu <- AddMetaData(object = A1_seu, metadata = combined_data, col.name = "category_to_colour_by")

####################################################################################################
## plot it 
SpatialPlot(object = A1_seu, group.by = "GOBP_LYMPHOCYTE_ACTIVATION", pt.size.factor = 2)

SpatialFeaturePlot(A1_seu, features = "GOBP_LYMPHOCYTE_ACTIVATION", pt.size.factor = 2.8)

###
# extract new seu metadata
meta_data_A1_pathways <- A1_seu@meta.data


################################################ same for B
## put it all together
###########################################################
## Seurat object
B1_seu <- readRDS("/Users/cathal.king/Documents/Projects/DL/Visium_Mouse/Seurat/Visium_mouse_Re_seurat_B1.rds")
SpatialPlot(object = B1_seu)
## pathway df
B_pathways <- read.csv(file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/GSVA_results/B1_Per_spot_activity_scores.csv", row.names = 1)
# colnames of the pathways df are the barcodes
# rownames of the pathways df are the pathways

# Remove the "-1" from the end of each column name
colnames(B1_seu) <- sub("-1$", "", colnames(B1_seu))

# test
# Remove the "_A1" from the end of each column name
colnames(B_pathways) <- sub("_B1$", "", colnames(B_pathways))
# Transpose and convert to a df
B_pathways_t <- as.data.frame(t(B_pathways))
# copy the rownames column (the 10x barcodes) to a different column in the pathways df
B_pathways_t <- B_pathways_t %>%
  mutate(Barcode = rownames(.))
# extract seu metadata
meta_data_B1 <- B1_seu@meta.data
# copy the rownames column (the 10x barcodes) to a different column in the meta_data
meta_data_B1_new <- meta_data_B1 %>%
  mutate(Barcode = rownames(.))

############# join
## now join the pathways df to the extracted meta-data df
meta_data_B1_joined <- left_join(x = meta_data_B1_new, y = B_pathways_t, by = "Barcode", keep=T)


### add that new meta data back into the Seurat object

## add the meta-data column with the matched barcodes to the seurat object. each pathway 1-by-1
B1_seu <- AddMetaData(object = B1_seu, metadata = meta_data_B1_joined, col.name = "GOBP_LYMPHOCYTE_ACTIVATION")

B1_seu <- AddMetaData(object = B1_seu, metadata = meta_data_B1_joined, col.name = "GOBP_B_CELL_ACTIVATION")
B1_seu <- AddMetaData(object = B1_seu, metadata = meta_data_B1_joined, col.name = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION")
B1_seu <- AddMetaData(object = B1_seu, metadata = meta_data_B1_joined, col.name = "GOBP_HORMONE_METABOLIC_PROCESS")

a <- SpatialFeaturePlot(B1_seu, features = "GOBP_LYMPHOCYTE_ACTIVATION", pt.size.factor = 2.8)
b <- SpatialFeaturePlot(B1_seu, features = "GOBP_B_CELL_ACTIVATION", pt.size.factor = 2.8)
c <- SpatialFeaturePlot(B1_seu, features = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION", pt.size.factor = 2.8)
d <- SpatialFeaturePlot(B1_seu, features = "GOBP_HORMONE_METABOLIC_PROCESS", pt.size.factor = 2.8)

##
# extract new seu metadata
meta_data_B1_pathways <- B1_seu@meta.data

write.csv(x = meta_data_B1_pathways, file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/Seurat_objects_with_pathways/meta_datas/B1_meta_pathway.csv")
###########################################################






################################################ same for C
## put it all together
###########################################################
## Seurat object
C1_seu <- readRDS("/Users/cathal.king/Documents/Projects/DL/Visium_Mouse/Seurat/Visium_mouse_Re_seurat_C1.rds")
SpatialPlot(object = C1_seu)
## pathway df
C_pathways <- read.csv(file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/GSVA_results/C1_Per_spot_activity_scores.csv", row.names = 1)
# colnames of the pathways df are the barcodes
# rownames of the pathways df are the pathways

# Remove the "-1" from the end of each column name
colnames(C1_seu) <- sub("-1$", "", colnames(C1_seu))

# test
# Remove the "_A1" from the end of each column name
colnames(C_pathways) <- sub("_C1$", "", colnames(C_pathways))
# Transpose and convert to a df
C_pathways_t <- as.data.frame(t(C_pathways))
# copy the rownames column (the 10x barcodes) to a different column in the pathways df
C_pathways_t <- C_pathways_t %>%
  mutate(Barcode = rownames(.))
# extract seu metadata
meta_data_C1 <- C1_seu@meta.data
# copy the rownames column (the 10x barcodes) to a different column in the meta_data
meta_data_C1_new <- meta_data_C1 %>%
  mutate(Barcode = rownames(.))

############# join
## now join the pathways df to the extracted meta-data df
meta_data_C1_joined <- left_join(x = meta_data_C1_new, y = C_pathways_t, by = "Barcode", keep=T)


### add that new meta data back into the Seurat object

## add the meta-data column with the matched barcodes to the seurat object. each pathway 1-by-1
C1_seu <- AddMetaData(object = C1_seu, metadata = meta_data_C1_joined, col.name = "GOBP_LYMPHOCYTE_ACTIVATION")

C1_seu <- AddMetaData(object = C1_seu, metadata = meta_data_C1_joined, col.name = "GOBP_B_CELL_ACTIVATION")
C1_seu <- AddMetaData(object = C1_seu, metadata = meta_data_C1_joined, col.name = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION")
C1_seu <- AddMetaData(object = C1_seu, metadata = meta_data_C1_joined, col.name = "GOBP_HORMONE_METABOLIC_PROCESS")

a <- SpatialFeaturePlot(C1_seu, features = "GOBP_LYMPHOCYTE_ACTIVATION", pt.size.factor = 2.8)
b <- SpatialFeaturePlot(C1_seu, features = "GOBP_B_CELL_ACTIVATION", pt.size.factor = 2.8)
c <- SpatialFeaturePlot(C1_seu, features = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION", pt.size.factor = 2.8)
d <- SpatialFeaturePlot(C1_seu, features = "GOBP_HORMONE_METABOLIC_PROCESS", pt.size.factor = 2.8)

##
# extract new seu metadata
meta_data_C1_pathways <- C1_seu@meta.data

write.csv(x = meta_data_C1_pathways, file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/Seurat_objects_with_pathways/meta_datas/C1_meta_pathway.csv")
###########################################################





################################################ same for D
## put it all together
###########################################################
## Seurat object
D1_seu <- readRDS("/Users/cathal.king/Documents/Projects/DL/Visium_Mouse/Seurat/Visium_mouse_Re_seurat_D1.rds")
SpatialPlot(object = D1_seu)
## pathway df
D_pathways <- read.csv(file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/GSVA_results/D1_Per_spot_activity_scores.csv", row.names = 1)
# colnames of the pathways df are the barcodes
# rownames of the pathways df are the pathways

# Remove the "-1" from the end of each column name
colnames(D1_seu) <- sub("-1$", "", colnames(D1_seu))

# test
# Remove the "_A1" from the end of each column name
colnames(D_pathways) <- sub("_D1$", "", colnames(D_pathways))
# Transpose and convert to a df
D_pathways_t <- as.data.frame(t(D_pathways))
# copy the rownames column (the 10x barcodes) to a different column in the pathways df
D_pathways_t <- D_pathways_t %>%
  mutate(Barcode = rownames(.))
# extract seu metadata
meta_data_D1 <- D1_seu@meta.data
# copy the rownames column (the 10x barcodes) to a different column in the meta_data
meta_data_D1_new <- meta_data_D1 %>%
  mutate(Barcode = rownames(.))

############# join
## now join the pathways df to the extracted meta-data df
meta_data_D1_joined <- left_join(x = meta_data_D1_new, y = D_pathways_t, by = "Barcode", keep=T)


### add that new meta data back into the Seurat object

## add the meta-data column with the matched barcodes to the seurat object. each pathway 1-by-1
D1_seu <- AddMetaData(object = D1_seu, metadata = meta_data_D1_joined, col.name = "GOBP_LYMPHOCYTE_ACTIVATION")

D1_seu <- AddMetaData(object = D1_seu, metadata = meta_data_D1_joined, col.name = "GOBP_B_CELL_ACTIVATION")
D1_seu <- AddMetaData(object = D1_seu, metadata = meta_data_D1_joined, col.name = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION")
D1_seu <- AddMetaData(object = D1_seu, metadata = meta_data_D1_joined, col.name = "GOBP_HORMONE_METABOLIC_PROCESS")

a <- SpatialFeaturePlot(D1_seu, features = "GOBP_LYMPHOCYTE_ACTIVATION", pt.size.factor = 2.8)
b <- SpatialFeaturePlot(D1_seu, features = "GOBP_B_CELL_ACTIVATION", pt.size.factor = 2.8)
c <- SpatialFeaturePlot(D1_seu, features = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION", pt.size.factor = 2.8)
d <- SpatialFeaturePlot(D1_seu, features = "GOBP_HORMONE_METABOLIC_PROCESS", pt.size.factor = 2.8)

cowplot::plot_grid(a,b,c,d, nrow = 2)
##
# extract new seu metadata
meta_data_D1_pathways <- D1_seu@meta.data

write.csv(x = meta_data_D1_pathways, file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/Seurat_objects_with_pathways/meta_datas/D1_meta_pathway.csv")
###########################################################




###### repeat this except for combined samples

################################################ same for C
## put it all together
###########################################################
## Seurat object
A1_seu <- readRDS("/Users/cathal.king/Documents/Projects/DL/Visium_Mouse/Seurat/Visium_mouse_Re_seurat_A1.rds")
SpatialPlot(object = A1_seu)
## pathway df
A_pathways <- read.csv(file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/GSVA_results/combined_analysis/A1_GSVA_activity_scores.csv", row.names = 1)
# colnames of the pathways df are the barcodes
# rownames of the pathways df are the pathways

# Remove the "-1" from the end of each column name
colnames(A1_seu) <- sub("-1$", "", colnames(A1_seu))

# test
# Remove the "_A1" from the end of each column name
colnames(A_pathways) <- sub("_A1$", "", colnames(A_pathways))
# Transpose and convert to a df
A_pathways_t <- as.data.frame(t(A_pathways))
# copy the rownames column (the 10x barcodes) to a different column in the pathways df
A_pathways_t <- A_pathways_t %>%
  mutate(Barcode = rownames(.))
# extract seu metadata
meta_data_A1 <- A1_seu@meta.data
# copy the rownames column (the 10x barcodes) to a different column in the meta_data
meta_data_A1_new <- meta_data_A1 %>%
  mutate(Barcode = rownames(.))

############# join
## now join the pathways df to the extracted meta-data df
meta_data_A1_joined <- left_join(x = meta_data_A1_new, y = A_pathways_t, by = "Barcode", keep=T)


### add that new meta data back into the Seurat object

## add the meta-data column with the matched barcodes to the seurat object. each pathway 1-by-1
A1_seu <- AddMetaData(object = A1_seu, metadata = meta_data_A1_joined, col.name = "GOBP_LYMPHOCYTE_ACTIVATION")
A1_seu <- AddMetaData(object = A1_seu, metadata = meta_data_A1_joined, col.name = "GOBP_B_CELL_ACTIVATION")
A1_seu <- AddMetaData(object = A1_seu, metadata = meta_data_A1_joined, col.name = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION")
A1_seu <- AddMetaData(object = A1_seu, metadata = meta_data_A1_joined, col.name = "GOBP_HORMONE_METABOLIC_PROCESS")

a <- SpatialFeaturePlot(A1_seu, features = "GOBP_LYMPHOCYTE_ACTIVATION", pt.size.factor = 2.8)
b <- SpatialFeaturePlot(A1_seu, features = "GOBP_B_CELL_ACTIVATION", pt.size.factor = 2.8)
c <- SpatialFeaturePlot(A1_seu, features = "GOBP_EPITHELIAL_CELL_DIFFERENTIATION", pt.size.factor = 2.8)
d <- SpatialFeaturePlot(A1_seu, features = "GOBP_HORMONE_METABOLIC_PROCESS", pt.size.factor = 2.8)

cowplot::plot_grid(a,b,c,d, nrow = 2)
##
# extract new seu metadata
meta_data_A1_pathways <- A1_seu@meta.data
# /Users/cathal.king/Documents/Projects/DL/featureCounts_clean/GSVA_results/combined_analysis
write.csv(x = meta_data_A1_pathways, file = "/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/GSVA_results/combined_analysis/seu_meta_datas/A1_meta_data.csv")
###########################################################




