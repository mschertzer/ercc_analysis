# Absolute RNA quantification using ERCC spike-ins

Megan D. Schertzer, [J. Mauro Calabrese](https://www.med.unc.edu/pharm/calabreselab/)

This ERCC analysis pipeline for absolute quantification of RNA-seq data is published alongside our [BioProtocols paper]().

## 1. RNA-Seq Processing (Fastq --> Bam)

### Getting Setup

**Analysis tools needed:**

I would recommend running these on a cluster

* STAR
* Samtools
* Subread- featureCounts

**Analysis ERCC files needed:** 

These files are located in the github repository or can be download from <https://www.thermofisher.com/order/catalog/product/4456739?SID=srch-srp-4456739#/4456739?SID=srch-srp-4456739> 

* ERCC92.fa
* ERCC92.gtf
* ERCC92_conc.txt

**Download genome annotation files for your cell type:**

You need a genome fasta file and a gtf file with gene features and coordinates. MM9 genome files were used in our analysis for the paper and are included in this github repository. Genome files compiled by the UCSC Genome Browser (Haeussler et al., 2019), which can be downloaded from Illumina’s iGenomes site: <https://support.illumina.com/sequencing/sequencing_software/igenome.html>.  

* mm9_genome.fa 
* mm9_genes.gtf

### *This next step is optional* 

If you want a dataset to test the pipeline on, you can download some of our data from GEO.

* Example file to download: VP160 noguide rep2 RNA-seq from GEO Series GSE118401
* Requires sratoolkit
* I would recommend running sratoolkit on a cluster

**Get file from GEO:**

```
fastq-dump SRR7685881

```

### Running the analysis pipeline

I have the file names that I used for my example file. When processing your own data, substitute the proper names and parameters. If you are on a cluster, these should be submitted as a job (sbatch, qsub, bsub, etc. depending on the platform). These should NOT be run on the login node.

**Build STAR genome that includes ERCC spike-in fasta and gtf:**

```
STAR 
  --runThreadN 8 
  --runMode genomeGenerate 
  --genomeDir mm9_STAR_ercc 
  --genomeFastaFiles mm9_genome.fa ERCC92.fa 
  --sjdbGTFfile mm9_genes.gtf ERCC92.gtf

```

**Align reads using STAR:**

```
STAR 
  --runThreadN 12 
  --genomeDir mm9_STAR_ercc 
  --readFilesIn SRR7685881.fastq 
  --outFileNamePrefix vprtta_rep2_
  
```

**Run samtools to filter for a mapping quality > 30:**

Filter, keep header information and output to a bam file.

```
samtools view -bhq 30 vprtta_rep2_Aligned.out.sam > vprtta_rep2_q30.bam

```

**Count reads in bam file:**

Get the read counts for the quality filtered alignemnt file. This step is important to calculate RPKM for the downstream analysis in Excel or R.

```
samtools view -c vprtta_rep2_q30.bam > vprtta_rep2_counts.txt
```

**Run featureCounts for ERCC spike-ins:**

```
featureCounts 
  -s 2
  -a ERCC92.gtf 
  -o vprtta_rep2_fc_ercc.txt 
  vprtta_rep2_q30.bam
  
```

**Run featureCounts for genes:**

```
featureCounts 
  -s 2 
  -a mm9_genes.gtf 
  -o vprtta_rep2_fc.txt
  vprtta_rep2_q30.bam
  
```
  
## 2. Analyze ERCC and gene featureCount files in either R

* If you would rather use excel, there is also an excel template in the github repository.
* If you are going to use the R script, keep reading. 

### Getting Setup

* I highly recommend using Rstudio to run the R script. It will be easier to navigate if the R script has issues installing and loading packages. 
* Using R studio, you can run R code line by line in the "Console" and you can run the entire R script on the command line in the "Terminal".
* If you have an Rstudio instance on your cluster, thats great! Use that.
* You can also download Rstudio locally on your computer.

Default parameters are based on what we used in the paper. Our trophoblast stem cells have about 30pg per cell (-p 30). We diluted the ERCC spike-in Mix 1 by 1:100 (-d 100). From this dilution, we added 2 microliters of ERCC mix (-e 2) to 1 microgram of total RNA from our cells (-r 1). 

### Running the R script

**Here is the generic command line prompt that you can use for your data. Just subtitute the files and parameters for your experiment:**

```
Rscript --vanilla ercc_analysis.R -p <RNA_PER_CELL_PICOGRAM> -d <DILUTION> -e <ERCC_MICROLITER> -r <RNA_MICROGRAM> <ERCC_FEATURECOUNTS_FILE> <GENE_FEATURECOUNTS_FILE> <BAM_READ_COUNT>

```

**If you want to test the R script on some test data, use the following command line argument to process provided featureCounts output using default parameters:**

These files can be found in the github repository.

```
Rscript --vanilla ercc_analysis.R -p 30 -d 100 -e 2 -r 1 vprtta_rna_ercc_fc.txt vprtta_rna_fc.txt 34820981

```

### *Help for running R script:*

There are **3 required inputs** for this R script (the order of the inputs matters!):

  1. ERCC featureCounts output (should only have 92 lines plus header information)
  2. Gene featureCounts output
  3. Read counts (This value can be optained using 'samtools view -c Aligned_reads.bam')
  
**Recommended parameters:**

-p RNA_PER_CELL_PICOGRAM, --rna_per_cell_picogram=RNA_PER_CELL_PICOGRAM

|       [Default 30] based on the cell type used in our paper. This is the amount of total RNA per cell for your cell type. This is cell type specific and must be calculated for your cell type as described in the paper.  
&nbsp;

**Additional parameters:**

Options:  

-d DILUTION, --dilution=DILUTION

|       The dilution of ercc spike-ins prepared prior to library prep. [Default 100] is     same as used in paper.  
&nbsp;

-e ERCC_MICROLITER, --ercc_microliter=ERCC_MICROLITER  

|       Microliter amount of ercc spike-in dilution added to rna before library prep.     [Default 2] is same as used in paper.  
&nbsp;

-r RNA_MICROGRAM, --rna_microgram=RNA_MICROGRAM  

|       Micrograms of starting total cellular RNA used in library prep. [Default 1] is same as used in paper.  
&nbsp;



  
