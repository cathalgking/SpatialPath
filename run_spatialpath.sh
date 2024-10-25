#!/bin/bash

spaceranger count --id=sample345 ./ \ #Output directory
                  --transcriptome=/home/jdoe/refdata/GRCh38-2020-A \ #Path to Reference
                  --probe-set=Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv \ #Path to probe set
                  --fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \ #Path to FASTQs
                  --sample=mysample \ #Sample name from FASTQ filename
                  --image=/home/jdoe/runs/images/sample345.tiff \ #Path to brightfield image
                  --slide=V19J01-123 \ #Slide ID
                  --area=A1 \ #Capture area
                  --localcores=8 \ #Allowed cores in localmode
                  --localmem=64  #Allowed memory (GB) in localmode

echo "SpaceRanger_Done"

cd /outs/ # Output directory of SpaceRanger

awk 'BEGIN{FS=","} NR>1 && $2==1 {print $1}' spatial/tissue_positions.csv > barcodes_in_tissue.csv

mkdir split

echo "start_samtools_split_process"

samtools view -u -D CB:barcodes_in_tissue.csv possorted_genome_bam.bam | samtools split -d CB -M 5000 --output-fmt bam -f 'split/%!.bam' -

echo "samtools_done_start_FeatureCounts"

featureCounts -T 28 -t exon -g gene_id -a  Mus_musculus.GRCm38.99.gtf -o featureCounts_SP.txt split/*bam

cut -f1-5 -d " " featureCounts_SP.txt | sed 's/split\///g' | cut -f 1,7- | grep -v ";" > filtered_feature_counts.txt

echo "start_Rscript_1"

Rscript path_to_SpatialPath/R_scripts/Spatial_pathR_step1_processing_lcpm.R

echo "start_Rscript_2"

Rscript path_to_SpatialPath/R_scripts/Spatial_pathR_step2_GSVA.R

cd ../

echo "start_Rscript_3"

Rscript /path_to_SpatialPath/R_scripts/Spatial_pathR_step3_Seurat.R

echo "Finished Running SpatialPath!"
