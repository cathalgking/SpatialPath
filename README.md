# SpatialPath

*Under construction...*

SpatialPath is a pipeline that provides pathway-level activity scores across spatial omics datasets. The goal of SpatialPath is to process your data so that you can get to the fun part of the analysis quicker.

The workflow for SpatialPath is shown below. Inputs include FASTQ's, along with histology images, which are put through numerous analyses including read mapping, gene expression quantification, spatial binning and pathway activity scoring. Outputs include spatial plots with pathway level scores overlayed across each sample. Data is conveniently formatted within commonly used object structures namely ```SpatialExperiment```, ```Seurat``` and ```Anndata``` which allows for further analyses or plot manipulation.


## Installation

To set up the **SpatialPath** environment and download the necessary reference files, run the installation script.

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html) installed  
- git


```bash
git clone https://github.com/cathalgking/SpatialPath.git
cd SpatialPath
./install_spatialpath.sh
```




## Usage

```bash
./spatialpath.sh -t 28 -i /path/to/spaceranger/output -s mouse
```

### Command-Line Options

- **`-t <threads>`**: Number of threads to use (default: 1)
- **`-i <input_folder>`**: Input folder (must be the output from Spaceranger)
- **`-s <species>`**: Species (`mouse` or `human`)
- **`-p <spatialpath_dir>`**: Location of SpatialPath installation
