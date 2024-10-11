
#dependancies edgeR, GSVA

raw_cnts = read.delim("filtered_feature_counts.txt",sep="\t",row.names = 1,check.names = FALSE)

## Gene ID conversion table
ens2gene = readRDS("/homes/feargal.ryan/programs/SpatialPath/Rdata/ens2gene.RDS")



### EdgeR ###
library(edgeR)
keep = apply(raw_cnts>10,1,sum) > 50
raw_cnts = raw_cnts[keep,]
y = DGEList(counts=raw_cnts)
y = calcNormFactors(y,method = "TMM")
lcpm <- cpm(y, log=TRUE,normalized.lib.sizes = TRUE)

## Relabel rows with gene name instead of ensembl gene ids
lcpm_counts_TMM_gene = as.data.frame(lcpm)
lcpm_counts_TMM_gene$Gene = ens2gene[rownames(lcpm),"external_gene_name"]
lcpm_counts_TMM_gene = lcpm_counts_TMM_gene[!is.na(lcpm_counts_TMM_gene$Gene),]
lcpm_counts_TMM_gene = lcpm_counts_TMM_gene[!duplicated(lcpm_counts_TMM_gene$Gene),]
rownames(lcpm_counts_TMM_gene) = lcpm_counts_TMM_gene$Gene
lcpm_counts_TMM_gene = lcpm_counts_TMM_gene[,colnames(lcpm_counts_TMM_gene) != "Gene"]
lcpm_counts_TMM_gene = as.matrix(lcpm_counts_TMM_gene)
#lcpm_counts_TMM_gene[is.na(lcpm_counts_TMM_gene)]
lcpm_counts_TMM_gene = t(t(lcpm_counts_TMM_gene))
saveRDS(lcpm_counts_TMM_gene,"counts4gsva.RDS")
