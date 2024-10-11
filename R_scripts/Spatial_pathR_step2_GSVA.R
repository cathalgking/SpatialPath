msigDB_C5 = readRDS("/homes/feargal.ryan/programs/SpatialPath/Rdata/msigDB_C5_mouse.RDS")
#This will need to be changed but for now it's fine. 
msigDB_C5 = msigDB_C5[c("GOBP_INNATE_IMMUNE_RESPONSE",
                        "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
                        "GOBP_HORMONE_METABOLIC_PROCESS",
                        "GOBP_MITOTIC_CELL_CYCLE",
                        "GOBP_AUTOPHAGIC_CELL_DEATH",
                        "GOBP_ATP_METABOLIC_PROCESS",
                        "GOBP_DEFENSE_RESPONSE_TO_TUMOR_CELL",
                        "GOBP_NEURONAL_SIGNAL_TRANSDUCTION",
                        "GOBP_CARDIOBLAST_DIFFERENTIATION",
                        "GOBP_PROTEIN_FOLDING")]

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

