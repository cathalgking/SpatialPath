# SpatialPath

##### A MetroMap of the SpatialPath pipeline

<img width="818" alt="SP_metroMap" src="https://github.com/cathalgking/SpatialPath/assets/32261323/04260943-b35f-46a6-accf-c97484dc621d">

### Quick start SpatialPath pipeline:

```ssh cathal.king@hpc-lin-cmp03```
```conda activate nf2```

![image](https://github.com/user-attachments/assets/422fa691-2cde-4404-8762-8c9bf7eb1bf3)


(singularity or a suitable conda env is usually required for running a nf-core pipeline)

Example:
```{bash}
nextflow run /homes/cathal.king/nf_pipelines/nf-core-spatialpath \
-profile singularity \
--reference /homes/cathal.king/References/refdata-gex-GRCh38-2020-A \
--probeset /homes/cathal.king/nf_pipelines/SP_yard_lung/visium_raw_data/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_probe_set.csv \
--input /homes/cathal.king/nf_pipelines/SP_yard_lung/lung_samplesheet.csv \
--outdir /homes/cathal.king/SP/testing/running_SP
```

Quick start notes:
* The samplesheet contains FASTQ's, image, sample etc.
* The 



Raw mouse CD40 data is here: ```/homes/feargal.ryan/data/spatial_path/SAGCQA0109_Visium```

```
All the files and everything for it can be found under /homes/feargal.ryan/programs/SpatialPath. There's a slurm script there called SpatialPath_v1.slurm
```

## **Aim**: Pathway Analysis for Spatial Transcriptomics data.

### The main steps are:
1. Align the data with SpaceRanger (***optional***)
2. Parse the output BAM file from SpaceRanger down to 1 BAM file per spatial spot. This can be done with samtools or the sinto package.
3. Next, featureCounts will be run on each per-spot BAM file to assign reads to Features.
4. Tidy the output from featureCounts so that each table has genes (ENSXXX) as the rownames and 10x barcodes as the colnames.
5. Use the GSVA R package to assign pathway scores per spot.
6. Merge GSVA output to spatial object and visualise results.

### Other potential names:
* SpatialGSVA
* PathMapR
* SpatialSignatureMapper

### Public Spatial Datasets

**10x Genomics**
https://www.10xgenomics.com/datasets

**10x Genomics Mouse Saggital Brain dataset as shown in the Seurat 10x Visium tutorial**
https://www.10xgenomics.com/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard

**MGI**
https://db.cngb.org/stomics/datasets/


## **Step by step notes on running the pipeline:**

### 1. Aligning the data with SpaceRanger
SpaceRanger is the 10x software used to align sequencing reads to Visium spots. It is a tool that is ran on the Linux/Unix command line and an example is shown below. The full list of options available for spaceranger count are given here --> https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/command-line-arguments

```{bash}
spaceranger count --id=V10U29_117_B1_Re_run \
 --transcriptome=/homes/cathal.king/yard/cell_r_metrics/refdata-gex-mm10-2020-A/ \
 --fastqs=/homes/john.salamon/visium/output/fastqc_output/outs/fastq_path/SAGCQA0109_Visium/21-01109/ \
 --sample=V10u29_B1 \
 --image=/data/dropbox/Visium/V10u29_Lynn_Group/V10u29_B1.tif \
 --slide=V10U29-117 \
 --area=B1 \
 --loupe-alignment=/homes/cathal.king/DL/Visium2/JSON_Alignments/V10U29-117-B1.json \
 --localcores=32 \
 --localmem=128
 --nosecondary
```

### 2. Parsing the BAM file

There is 2 methods that achieve the exact same thing for this step. 
The files that are required are the BAM file from SpaceRanger which is found in the ```/outs/``` folder. 

Also the spatial barcode file which is a text file of spatial barcodes for each Visium spot. This can be gotten from the first column of the file ```/outs/spatial/tissue_positions.csv``` from SpaceRanger outputs. An example of this file is found here --> ```/homes/cathal.king/LB/AR_kd/SpatialPath/spbar_lists/siCon3_M14_C1_spbars_new.txt```

The ```tissue_positions.csv``` file contains 6 columns but only the first 2 are needed for SpatialPath. The first column contains the 10x spot barcodes and the 2nd column (in_tissue) contains 1/0 based on whether the Visium spot is below tissue. Only the spots associated with tissue are required for SpatialPath (most uses) so subset down that file to extract only the spot barcodes with a "1" entry:
```awk 'BEGIN{FS=","} NR>1 && $2==1 {print $1}' tissue_positions.csv > barcodes_in_tissue.csv```

** The ```tissue_positions.csv``` file is sometimes named as ```tissue_positions_list.csv``` depending on the version on SpaceRanger being used.

For sinto, a file or comma-separated list of cell barcodes will suffice. However, a second column is required which is supposed to be the group that each cell belongs to. For this column, I just copied the 2nd column. This might be something else to optimize.


#### Method 1: samtools 
The version of samtools used (**v1.19**) is important here as the only versions dont have that functionality.
One of the outputs of SpaceRanger is the BAM file which contains the aligned reads. This file can be parsed down to only include reads from a single spatial spot or multiple spots. Each entry of the output BAM file contains a ```CB``` flag and this is used for this process. The user only needs to provide a list of spatial barcodes (spbars.txt) in a text file and the ```samtools view``` command is then as follows:



```{bash}
for i in $(cat PATH/spbars.txt) ; do samtools view -b -d CB:$i /PATH/outs/possorted_genome_bam.bam > $i.bam ; done
```

* ```-b``` flag specifies the output as BAM
* ```-d``` flag says to only output alignments that match from the ```CB``` entries

Feargal suggested this process could be sped up with the following:

```{bash}
cat PATH/spbars.txt | xargs -I {} -P <number_of_parallel_processes> sh -c 'samtools view -b -d CB:$1 /PATH/outs/possorted_genome_bam.bam > $1.bam' -- {}
```


#### Method 2: sinto

It turns out that samtools takes a long time to Parse each BAM file (even on the SAHMRI HPC) and is memory intensive. This is because (I think) it is reading down through each BAM file read-by-read and that can be a slow process.

An alternative package is ```sinto``` which essentially does the same thing but much quicker. Sinto is also ran on the command line and the function to parse down a BAM file is ```filterbarcodes```. An example command would look like:

First, might have to copy the barcodes file to have duplicate columns with:

```awk '{print $0,$NF}' barcodes_in_tissue.csv > new_barcodes.csv```


```{bash}
  sinto filterbarcodes \
  --bam /PATH/possorted_genome_bam.bam \
  --nproc 10 \
  --cells /homes/cathal.king/LB/AR_kd/SpatialPath/spbar_lists/siCon1_F28_C1_spbars_new.txt \
  --outdir /homes/cathal.king/LB/AR_kd/SpatialPath/sinto_all/siCon1_F28_C1_spbars
```
sinto is on the SAHMRI HPC and can be accessed by

```{bash}
conda activate /cancer/storage/tools/conda_apps/sinto/env/
#
sinto -v # check installation
```

sinto --> https://github.com/timoast/sinto

#### Method 3: script that Paul wrote

```/cancer/storage/SAGC/scratch/split_bam_by_tag_noSort.py```


### 3. Assign reads to Features with featureCounts

Use featureCounts to assign reads to features. This command is ran on each BAM file (1 per-spot) and the output will be a txt file.

```{bash}
featureCounts -T 32 -t exon -g gene_id 
-a ~/databases/grcm38/Mus_musculus.GRCm38.99.gtf 
-o featureCounts_output.txt /homes/PATH/to/Parsed_BAMS/A1_spots/*/*bam
```

### 4. Clean the output from featureCounts

This step might vary depending on the data but at the end the featureCounts table should have genes (ENSXXX) as the rownames and 10x spatial barcodes (e.g AACACTTGGCAAGGAA-1) as the colnames.

```{bash}
# then remove the 2nd to the 6th column which has data that is not needed
cut -f 1,7- output.txt > output.txt

# remove the first row from the file which just has info on the fCount command that was used
tail -n +2 featureCounts_output.txt > output.txt

# remove the PATH string from the beginning of each column name. 
#sed 's#/homes/cathal.king/LB/AR_kd/SpatialPath/spbars/siAR3_M14_D1_spbars/bams_perSpot/##g' output.txt > output.txt
# and also remove the ".bam" from each colname
# the below sed line will do it all in 1 line
sed 's#/homes/cathal.king/LB/AR_kd/SpatialPath/spbars/siAR3_M14_D1_spbars/bams_perSpot/##g; s/\.bam//g' output.txt > featureCounts_Cleaned.txt
```

### 5. Use the GSVA R package to assign pathway scores per spot.

Import the featureCount cleaned output into R. Run a GSVA analysis on that data as outlined in the script ```GSVA.R```.

### 6. Merge GSVA output to spatial object and visualise results.

After running the GSVA R script, there will be a csv file containing an activity pathway score for each Visium spot. An example of this data table looks like the below. This table might have been cleaned up a little for example with the names of the pathways.

<img width="533" alt="Screenshot 2024-01-18 at 10 34 02 am" src="https://github.com/cathalgking/SpatialPath/assets/32261323/eb0358bd-b141-43e7-adc4-051792ad996d">



This data is then added to the meta-data slot of the Spatial object containing the Visium sample. It is essentially joining two columns together based on a variable which is the 10x barcode column (colnames) in the Spatial R object.



* Ensure that the GSVA scores are joined to the object wrt the 10x barcodes (colnames). It is **vital** that this step is done correctly because otherwise, GSVA scores will be getting assigned to the wrong Visium spots!
* The Spatial R object could be a Seurat object or a SpatialExperiment (SPE) object. The principle is the same for both but the methods differ slightly. Let us work with a Seurat object now for simplicity.

***Example meta-data of a Seurat object***

<img width="519" alt="Screenshot 2024-01-18 at 10 38 03 am" src="https://github.com/cathalgking/SpatialPath/assets/32261323/6ae038df-a014-4ad5-b0c1-a2945b403148">


Steps to add GSVA scores

1-Read in Visium data and create a Seurat object

2-Extract the meta-data in the Seu object to a seperate data-frame

3-Read in the GSVA scores df

4-Append the GSVA scores df to the extracted meta-data df.

5-Add that new df back into the Seurat object meta-data.

```{R}
library(Seurat)
library(SeuratObject)

# the directory containing the SpaceRanger "/outs/" folder
data_dir <- "/Users/cathal.king/Documents/Projects/LB/AR_kd/new_outs/new_outs/siAR1_F28_D1/"

# create a Seurat object
seu_object <- Seurat::Load10X_Spatial(data.dir = data_dir)

# ensure Seurat object is loaded correctly. It should say something like "An object of class Seurat".
seu_object

## Examine the meta-data of this Seurat object. This can be accessed multiple ways. 
# start with
head(seu_object) # this should look like the above screen-shot

# write to a seperate
seu_meta <- seu_object@meta.data

# read in GSVA table
GSVA_paths <- read.csv("/Users/cathal.king/Documents/Projects/LB/AR_kd/pathways/GSVA_scores_per_spot_cleaned/siAR1_F28_D1_gsva_transposed.csv")
GSVA_paths <- GSVA_paths[,1:2] # only keep 1 pathway for simplicity

# transpose pathways df
#GSVA_paths_t <- t(GSVA_paths)

# copy the rownames column (the 10x barcodes) to a different column in the meta_data. This is so we can check that the join was done correctly.
seu_meta_new <- seu_meta %>%
  mutate(Barcode = rownames(.))

## now join the pathways df to the extracted meta-data df
meta_data_joined <- left_join(x = seu_meta_new, y = GSVA_paths, by = "Barcode", keep=T)

### add that new meta data back into the Seurat object

## add the meta-data column with the matched barcodes to the seurat object. each pathway 1-by-1
seu_object <- AddMetaData(object = seu_object, metadata = meta_data_joined, col.name = "GOBP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY")

# check the meta-data slot again. It should now contain the pathway column.
head(seu_object)
```

* This step can probably optimised a lot more. Just exercising caution for now in going step-by-step.
* Sometimes the order of the barcodes is the same from the Seu object and the GSVA object and a column bind would do the same thing.
* The ```AddMetaData``` function is basically a column bind and it does not re-order the column to match when joining.

### 7. Plot the result
* The pathway scores should now be added to the spatial object.
* Plot the result by referring to the ```group.by``` argument in the ```SpatialPlot()``` function.
* This is a Seurat function.

```{r}
## plot it 
SpatialPlot(object = A1_seu, group.by = "GOBP_LYMPHOCYTE_ACTIVATION", pt.size.factor = 2)

SpatialFeaturePlot(A1_seu, features = "GOBP_LYMPHOCYTE_ACTIVATION", pt.size.factor = 2.8)
```



