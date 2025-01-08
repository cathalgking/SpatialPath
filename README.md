# SpatialPath

SpatialPath is a pipeline that provides pathway-level activity scores across spatial omics datasets. 

The workflow for SpatialPath is shown below. Inputs include FASTQ's, along with histology images, which are put through numerous analyses including read mapping, gene expression quantification, spatial binning and pathway activity scoring. Outputs include spatial plots with pathway level scores overlayed across each sample.


## Installation Instructions

To set up the **SpatialPath** environment and download the necessary reference files, run the installation script.

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html) installed  
- git


```bash
git clone https://github.com/cathalgking/SpatialPath.git
cd SpatialPath
./install_SpatialPath.sh
```




## Usage

```bash
./spatialpath.sh -t 28 -i /path/to/spaceranger/output -s mouse -p /path/to/spatialpath
```

### Command-Line Options

- **`-t <threads>`**: Number of threads to use (default: 1)
- **`-i <input_folder>`**: Input folder (must be the output from Spaceranger)
- **`-s <species>`**: Species (`mouse` or `human`)
- **`-p <spatialpath_dir>`**: Location of SpatialPath installation


## Example

### Dataset
```SpatialPath``` is demonstrated below on a 10x Visium CytAssist experiment<sup>1</sup> on fresh frozen mouse brain tissue. The tissue section contained 4,298 detected spatial spots (voxels), with a median of 19,627 UMI counts and 6,178 genes per spot. Sequencing was conducted on an Illumina NovaSeq, achieving a sequencing depth of 171,410,389 reads with 38.4% saturation. Multiple pathways of interest can be input to ```SpatialPath```. For this dataset, all GO pathways were input including common neurological pathways such as Neurotransmitter secretion.

![Alt text](svg/mouse_brain_image.svg)


### Outputs

```SpatialPath``` has multiple outputs which allow for easy visualisation and customisation. One spatial plot visualisation per input pathway will be outputted in PDF format. Raw data is stored within commonly used data objects namely ```SpatialExperiment```, ```Seurat``` and ```Anndata``` for easy access and customisation.  Below are four example SpatialPath plots, displaying one GO pathway per plot, for the example mouse brain dataset.

![Alt text](svg/mouse_brain_NeuroinflammatoryResponse.svg)

![Alt text](svg/mouse_brain_NeurotransmitterSecretion.svg)

![Alt text](svg/mouse_brain_LIPID_METABOLIC_PROCESS.svg)

![Alt text](svg/mouse_brain_OXIDATIVE_PHOSPHORYLATION.svg)


## References

1 - https://www.10xgenomics.com/datasets/fresh-frozen-visium-on-cytassist-mouse-brain-probe-based-whole-transcriptome-profiling-2-standard
