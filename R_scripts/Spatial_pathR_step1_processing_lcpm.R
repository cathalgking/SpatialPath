
#dependancies edgeR, GSVA

raw_cnts = read.delim("filtered_feature_counts.txt",sep="\t",row.names = 1,check.names = FALSE)


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
path_to_installation <- args[1]
species <- args[2]

# Construct the species-specific path
path_to_rdata <- file.path(path_to_installation, "R_scripts", species)

ens2gene <- readRDS(file.path(path_to_rdata, "Rdata", "ens2gene.RDS"))


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
