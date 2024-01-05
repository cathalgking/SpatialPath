# SpatialPath

* Aim: to visualise pathway level activity scores across a spatial transcriptomics image.

## Other potential names:
* SpatialGSVA
* PathMapR
* SpatialSignatureMapper

Colnames in a Seurat object are the 10x spatial barcodes for example TTGACAGGAGCTCCCG-1

### Extract colnames from a Seurat object
# just 1 Seurat object
spbar_siar3 <- colnames(siAR3_M14_D1)
# write to a file
write.table(spbar_siar3, file = "/PATH/SpatialPath/test_spbar7.txt", row.names = F, quote = F, col.names = F)
