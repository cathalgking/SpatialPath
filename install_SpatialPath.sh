#!/bin/bash

# Set environment name and GTF directory
ENV_NAME="SpatialPath_v0.1"
GTF_DIR="references"

# URLs for GTF files
MOUSE_GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz"
HUMAN_GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"

# Usage function to display help
usage() {
  echo "Usage: $0 [-h]"
  echo ""
  echo "Options:"
  echo "  -h  Display this help message"
  exit 0
}

# Parse command-line options
while getopts "h" opt; do
  case ${opt} in
    h) usage ;;
    *) usage ;;
  esac
done

# Create the conda environment
echo "Creating conda environment: $ENV_NAME"
conda env create -n "$ENV_NAME" -f SpatialPath.yml || { echo "Failed to create conda environment"; exit 1; }

# Activate the environment
echo "Activating conda environment: $ENV_NAME"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME" || { echo "Failed to activate conda environment"; exit 1; }

# Create the GTF directory if it doesn't exist
mkdir -p "$GTF_DIR"

# Download GTF files
echo "Downloading GTF files..."
curl -L "$MOUSE_GTF_URL" -o "$GTF_DIR/Mus_musculus.GRCm39.113.gtf.gz" || { echo "Failed to download mouse GTF"; exit 1; }
curl -L "$HUMAN_GTF_URL" -o "$GTF_DIR/Homo_sapiens.GRCh38.113.gtf.gz" || { echo "Failed to download human GTF"; exit 1; }

# Decompress the GTF files
echo "Decompressing GTF files..."
gunzip -f "$GTF_DIR/Mus_musculus.GRCm39.113.gtf.gz" || { echo "Failed to decompress mouse GTF"; exit 1; }
gunzip -f "$GTF_DIR/Homo_sapiens.GRCh38.113.gtf.gz" || { echo "Failed to decompress human GTF"; exit 1; }

