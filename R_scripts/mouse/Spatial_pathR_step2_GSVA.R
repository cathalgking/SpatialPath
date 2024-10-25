# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the path is provided, else use a default value
if (length(args) == 0) {
  path_to_installation <- "Rdata"
} else {
  path_to_installation <- args[1]
}

path_to_installation = paste(path_to_installation,"/R_scripts/mouse/"


# Read in the MSigDB C5 data using the provided installation path
msigDB_C5 <- readRDS(file.path(path_to_installation, "msigDB_C5_mouse.RDS"))

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

