
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
path_to_installation <- args[1]
species <- args[2]

# Construct the species-specific path
path_to_rdata <- file.path(path_to_installation, "R_scripts", species)

ens2gene <- readRDS(file.path(path_to_rdata, "Rdata", "ens2gene.RDS"))
msigDB_C5 <- readRDS(file.path(path_to_rdata, "Rdata", "mSigDB.RDS"))

# Read in the pathways to plot
paths_to_plot <- read.csv(file.path(path_to_installation, "GO_terms.csv"), header = TRUE)




paths_to_plot = paths_to_plot[paths_to_plot$Plot==TRUE,1]

paths_to_plot = paths_to_plot[paths_to_plot %in% names(msigDB_C5)]
#This will need to be changed but for now it's fine. 
msigDB_C5 = msigDB_C5[paths_to_plot]

my_modules = list()

for (module in names(msigDB_C5)){
my_modules[[module]] = msigDB_C5[[module]]
}


lcpm_counts_TMM_gene = readRDS("counts4gsva.RDS")

### GSVA to calculate module activity scores
library(GSVA)
all_btm.gsva = gsva(expr=lcpm_counts_TMM_gene, gset.idx.list = (my_modules),min.sz=4)
per_spot_scores = t(all_btm.gsva)
per_spot_scores = data.frame(Barcode=rownames(per_spot_scores),per_spot_scores)
rownames(per_spot_scores) = NULL
saveRDS(per_spot_scores,"per_spot_scores.RDS")

