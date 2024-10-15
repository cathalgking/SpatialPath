# SpatialPath

*Under construction...*

SpatialPath is a pipeline that provides pathway-level activity scores across spatial omics datasets. The goal of SpatialPath is to process your data so that you can get to the fun part of the analysis quicker.

The workflow for SpatialPath is shown below. Inputs include FASTQ's, along with histology images, which are put through numerous analyses including read mapping, gene expression quantification, spatial binning and pathway activity scoring. Outputs include spatial plots with pathway level scores overlayed across each sample. Data is conveniently formatted within commonly used object structures namely ```SpatialExperiment```, ```Seurat``` and ```Anndata```.



## Dependencies

- [SpaceRanger](https://www.10xgenomics.com/support/software/space-ranger/downloads/space-ranger-installation)
- Seurat
- R
- GSVA

## Installation

For installation of SpaceRanger, please refer to the [10X Genomics installation guide](https://www.10xgenomics.com/support/software/space-ranger/downloads/space-ranger-installation).

To install the remaining dependencies, you can create a conda environment using the provided environment file:

```bash
conda env create -f SpatialPath.yml
```
