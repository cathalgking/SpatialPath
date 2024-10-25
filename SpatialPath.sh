#!/bin/bash

# Default values
threads=1
input_folder=""
species=""
spatialpath_dir=""

# Usage function to display help message
usage() {
  echo "Usage: $0 -t <threads> -i <input_folder> -s <species> -p <spatialpath_directory> [-h]"
  echo ""
  echo "Options:"
  echo "  -t  Number of threads (default: 1)"
  echo "  -i  Input folder (required, output of spaceranger)"
  echo "  -s  Species (mouse or human, required)"
  echo "  -p  Path to spatialpath installation directory (required)"
  echo "  -h  Display this help message"
  echo ""
  exit 0
}

# Parse command-line options
while getopts "t:i:s:p:h" opt; do
  case ${opt} in
    t) threads=${OPTARG} ;;
    i) input_folder=${OPTARG} ;;
    s) species=${OPTARG} ;;
    p) spatialpath_dir=${OPTARG} ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Check if required options are provided
if [[ -z "$input_folder" || -z "$species" || -z "$spatialpath_dir" ]]; then
  echo "Error: Input folder, species, and spatialpath directory are required."
  usage
fi

# Validate species input
if [[ "$species" != "mouse" && "$species" != "human" ]]; then
  echo "Error: Species must be either 'mouse' or 'human'."
  usage
fi

# Select GTF file based on species
gtf_file=""
if [[ "$species" == "mouse" ]]; then
  gtf_file="$spatialpath_dir/references/Mus_musculus.GRCm38.99.gtf"
elif [[ "$species" == "human" ]]; then
  gtf_file="$spatialpath_dir/references/Homo_sapiens.GRCh38.99.gtf"
fi

# Print settings for confirmation
echo "Running SpatialPath with the following parameters:"
echo "  Threads: $threads"
echo "  Input Folder: $input_folder"
echo "  Species: $species"
echo "  GTF File: $gtf_file"
echo "  SpatialPath Directory: $spatialpath_dir"

# Generate barcodes_in_tissue.csv
awk 'BEGIN{FS=","} NR>1 && $2==1 {print $1}' "$input_folder/spatial/tissue_positions.csv" > barcodes_in_tissue.csv

# Create split directory
mkdir -p split

# Run samtools split
samtools view -u -D CB:barcodes_in_tissue.csv possorted_genome_bam.bam | \
samtools split -d CB -M 5000 --output-fmt bam -f 'split/%!.bam' -

# Run featureCounts
featureCounts -T "$threads" -t exon -g gene_id -a "$gtf_file" -o featureCounts_SP.txt split/*bam

# Filter featureCounts output
cut -f1-5 -d " " featureCounts_SP.txt | sed 's/split\///g' | cut -f 1,7- | grep -v ";" > filtered_feature_counts.txt

# Run R scripts
Rscript "$spatialpath_dir/R_scripts/Spatial_pathR_step1_processing_lcpm.R"
Rscript "$spatialpath_dir/R_scripts/Spatial_pathR_step2_GSVA.R"

cd ../
Rscript "$spatialpath_dir/R_scripts/Spatial_pathR_step3_Seurat.R"
