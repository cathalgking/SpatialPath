## GSVA pathway activity scores

### for Visium mouse data
library(readxl)
library(ggplot2)
library(gridExtra)
library(reshape)

#meta_data = as.data.frame(read_xlsx("Meta data.xlsx",sheet = "Sheet1",))
#rownames(meta_data) = meta_data$SampleID
#/Users/cathal.king/Documents/Projects/DL/featureCounts_clean
#raw_cnts = read.delim("/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/B1_featureCounts_clean.txt",sep="\t",row.names = 1,check.names = FALSE)
A1_raw_cnts = read.delim("/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/A1_featureCounts_clean.txt",sep="\t",row.names = 1,check.names = FALSE)
B1_raw_cnts = read.delim("/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/B1_featureCounts_clean.txt",sep="\t",row.names = 1,check.names = FALSE)
C1_raw_cnts = read.delim("/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/C1_featureCounts_clean.txt",sep="\t",row.names = 1,check.names = FALSE)
D1_raw_cnts = read.delim("/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/D1_featureCounts_clean.txt",sep="\t",row.names = 1,check.names = FALSE)

raw_cnts = cbind(A1_raw_cnts,B1_raw_cnts,C1_raw_cnts,D1_raw_cnts)
#meta_data = meta_data[rownames(meta_data) %in% colnames(raw_cnts),]
#raw_cnts = raw_cnts[,rownames(meta_data)]
#colnames(raw_cnts) == rownames(meta_data)


#attributes = listAttributes(ensembl)
ensids=rownames(raw_cnts)
table(duplicated(ensids))

## BiomaRt ##
library(biomaRt)
ensembl=useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
ens2gene = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype','description','chromosome_name',"entrezgene_id","description"), 
                 filters = 'ensembl_gene_id', 
                 values = ensids, 
                 mart = ensembl)
ens2gene = ens2gene[!duplicated(ens2gene$ensembl_gene_id),]
rownames(ens2gene) = ens2gene$ensembl_gene_id


### EdgeR ###
library(edgeR)
keep = apply(raw_cnts>10,1,sum) > 50
table(keep)
raw_cnts = raw_cnts[keep,]
y = DGEList(counts=raw_cnts)
y = calcNormFactors(y,method = "TMM")
lcpm <- cpm(y, log=TRUE,normalized.lib.sizes = TRUE)


## C5 ontology Gene Sets
library(msigdbr)
msig.df = as.data.frame(msigdbr(species = "Mus musculus",category = "C5"))
msigDB = by(msig.df$gene_symbol,
            msig.df$gs_name,
            function(x) as.character(x))

msigDB = msigDB[c("GOBP_B_CELL_ACTIVATION","GOBP_EPITHELIAL_CELL_DIFFERENTIATION","GOBP_HORMONE_METABOLIC_PROCESS","GOBP_LYMPHOCYTE_ACTIVATION")]
names(msigDB)  

class(msigDB)

msigDB = as.list(msigDB)
class(msigDB)

my_modules = list()
my_modules[["GOBP_LYMPHOCYTE_ACTIVATION"]] = msigDB$GOBP_LYMPHOCYTE_ACTIVATION
my_modules[["GOBP_B_CELL_ACTIVATION"]] = msigDB$GOBP_B_CELL_ACTIVATION
my_modules[["GOBP_EPITHELIAL_CELL_DIFFERENTIATION"]] = msigDB$GOBP_EPITHELIAL_CELL_DIFFERENTIATION
my_modules[["GOBP_HORMONE_METABOLIC_PROCESS"]] = msigDB$GOBP_HORMONE_METABOLIC_PROCESS

## Relabel rows with gene name instead of ensembl gene ids
lcpm_counts_TMM_gene = as.data.frame(lcpm)
lcpm_counts_TMM_gene$Gene = ens2gene[rownames(lcpm),"external_gene_name"]
lcpm_counts_TMM_gene = lcpm_counts_TMM_gene[!is.na(lcpm_counts_TMM_gene$Gene),]
lcpm_counts_TMM_gene = lcpm_counts_TMM_gene[!duplicated(lcpm_counts_TMM_gene$Gene),]
rownames(lcpm_counts_TMM_gene) = lcpm_counts_TMM_gene$Gene
lcpm_counts_TMM_gene = lcpm_counts_TMM_gene[,colnames(lcpm_counts_TMM_gene) != "Gene"]
lcpm_counts_TMM_gene = as.matrix(lcpm_counts_TMM_gene)
lcpm_counts_TMM_gene[is.na(lcpm_counts_TMM_gene)]
lcpm_counts_TMM_gene = t(t(lcpm_counts_TMM_gene))
lcpm_counts_TMM_gene = lcpm_counts_TMM_gene[rownames(lcpm_counts_TMM_gene) %in% unlist(my_modules),]

library(GSVA)
### GSVA to calculate module activity scores
all_btm.gsva = gsva(expr=lcpm_counts_TMM_gene, gset.idx.list = (my_modules),min.sz=4)
colnames(all_btm.gsva)

A1_gsva = all_btm.gsva[,grepl("_A1",colnames(all_btm.gsva))]
B1_gsva = all_btm.gsva[,grepl("_B1",colnames(all_btm.gsva))]
C1_gsva = all_btm.gsva[,grepl("_C1",colnames(all_btm.gsva))]
D1_gsva = all_btm.gsva[,grepl("_D1",colnames(all_btm.gsva))]
#
write.csv(D1_gsva, file="/Users/cathal.king/Documents/Projects/DL/featureCounts_clean/GSVA_results/combined_analysis/D1_GSVA_activity_scores.csv")

