## Test SpatialPath (after SpaceRanger) on a larger dataset
* The mouse dataset from the Seurat tutorial

BAM file (22G) is here --> ```/homes/cathal.king/SP_yard2/test_data```

Compare 2 methods to Parse BAM file

1. Original samtools for loop ```for i in $(cat PATH/spbars.txt) ; do samtools view -b -d CB:$i /PATH/outs/possorted_genome_bam.bam > $i.bam ; done```
2. New Feargal method ```cat PATH/spbars.txt | xargs -I {} -P <number_of_parallel_processes> sh -c 'samtools view -b -d CB:$1 /PATH/outs/possorted_genome_bam.bam > $1.bam' -- {}```
