# Script: gene_lookup_tables.R
# Author: Feargal Ryan
# GitHub: https://github.com/feargalr
# 
# Description:
# This R script uses the biomaRt package to create lookup tables for human and mouse genes.
# It connects to the Ensembl database to retrieve gene IDs and additional annotations such as 
# gene names, biotypes, descriptions, and Entrez IDs. The resulting tables are saved as RDS files 
# for future use. This script is useful for converting between gene ID types across species.

## BiomaRt ##
library(biomaRt)

# ---------------------------------------
# Creating human gene lookup table
# ---------------------------------------
# Connect to the Ensembl database for human genes
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get Ensembl gene IDs
human_gene_ids <- getBM(attributes = c("ensembl_gene_id"), mart = ensembl_human)

# Create a lookup table with various gene annotations
human_ens2gene <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 
                 'description', 'chromosome_name', 'entrezgene_id'), 
  filters = 'ensembl_gene_id', 
  values = human_gene_ids, 
  mart = ensembl_human
)

# Remove duplicate entries and set row names
human_ens2gene <- human_ens2gene[!duplicated(human_ens2gene$ensembl_gene_id), ]
rownames(human_ens2gene) <- human_ens2gene$ensembl_gene_id

# Save the human lookup table
saveRDS(human_ens2gene, "human_ens2gene.RDS")

# ---------------------------------------
# Creating mouse gene lookup table
# ---------------------------------------
# Connect to the Ensembl database for mouse genes
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get Ensembl gene IDs
mouse_gene_ids <- getBM(attributes = c("ensembl_gene_id"), mart = ensembl_mouse)

# Create a lookup table with various gene annotations
mouse_ens2gene <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 
                 'description', 'chromosome_name', 'entrezgene_id'), 
  filters = 'ensembl_gene_id', 
  values = mouse_gene_ids, 
  mart = ensembl_mouse
)

# Remove duplicate entries and set row names
mouse_ens2gene <- mouse_ens2gene[!duplicated(mouse_ens2gene$ensembl_gene_id), ]
rownames(mouse_ens2gene) <- mouse_ens2gene$ensembl_gene_id

# Save the mouse lookup table
saveRDS(mouse_ens2gene, "mouse_ens2gene.RDS")