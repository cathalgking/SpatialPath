# SpatialPath

## **Aim**: to visualise pathway level activity scores across a spatial transcriptomics image.

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

### Public Dataset

**10x Genomics**
https://www.10xgenomics.com/datasets

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
```

### 2. Parsing the BAM file

There is 2 methods that achieve the exact same thing for this step. 
The files that are required are the BAM file from SpaceRanger which is found in the ```/outs/``` folder. 

Also the spatial barcode file which is a text file of spatial barcodes for each Visium spot. This can be gotten from the first column of the file ```/outs/spatial/tissue_positions.csv``` from SpaceRanger outputs. An example of this file is found here --> ```/homes/cathal.king/LB/AR_kd/SpatialPath/spbar_lists/siCon3_M14_C1_spbars_new.txt```

For sinto, a file or comma-separated list of cell barcodes will suffice. However, a second column is required which is supposed to be the group that each cell belongs to. For this column, I just copied the 2nd column. This might be something else to optimize.


#### Method 1: samtools 
One of the outputs of SpaceRanger is the BAM file which contains the aligned reads. This file can be parsed down to only include reads from a single spatial spot or multiple spots. Each entry of the output BAM file contains a ```CB``` flag and this is used for this process. The user only needs to provide a list of spatial barcodes (spbars.txt) in a text file and the ```samtools view``` command is then as follows:

```{bash}
for i in $(cat PATH/spbars.txt) ; do samtools view -b -d CB:$i /PATH/outs/possorted_genome_bam.bam > $i.bam ; done
```

* ```-b``` flag specifies the output as BAM
* ```-d``` flag says to only output alignments that match from the ```CB``` entries

#### Method 2: sinto

It turns out that samtools takes a long time to Parse each BAM file (even on the SAHMRI HPC) and is memory intensive. This is because (I think) it is reading down through each BAM file read-by-read and that can be a slow process.

An alternative package is ```sinto``` which essentially does the same thing but much quicker. Sinto is also ran on the command line and the function to parse down a BAM file is ```filterbarcodes```. An example command would look like:

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
# remove the first row from the file which just has info on the fCount command that was used
tail -n +2 featureCounts_output.txt > output.txt

# then remove the 2nd to the 6th column which has data that is not needed
cut -f 1,7- output.txt > output.txt

# remove the PATH string from the beginning of each column name. 
#sed 's#/homes/cathal.king/LB/AR_kd/SpatialPath/spbars/siAR3_M14_D1_spbars/bams_perSpot/##g' output.txt > output.txt
# and also remove the ".bam" from each colname
sed 's#/homes/cathal.king/LB/AR_kd/SpatialPath/spbars/siAR3_M14_D1_spbars/bams_perSpot/##g; s/\.bam//g' output.txt > featureCounts_Cleaned.txt
```

### 5. Use the GSVA R package to assign pathway scores per spot.

Import the featureCount cleaned output into R. Run a GSVA analysis on that data as outlined in the script ```GSVA.R```.

### 6. 

At the end of the GSVA R script, there will be a csv file which has an activity score (based on the chosen pathway) for each Visium spot. This file is then added to the Seurat object (or Spatial object) making sure that the barcodes are matched up correctly.

```{R}
# read in SpaceRanger output and construct a Seurat object
data <- 
```

