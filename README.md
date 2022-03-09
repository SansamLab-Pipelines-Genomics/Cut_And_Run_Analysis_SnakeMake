# Cut&Run and ChipSeq Pipeline


## Project Description
This is a workflow for processing Cut and Run (modified ChiP-seq) data. Original workflow courtesy of the Sansam Lab.

This README details the individual steps of the pipeline. To run the automated pipeline see the README about the bash script.

## Table of Contents

* [Requirements](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#requirements)
* [Description of individual steps in pipeline](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#run-the-pipeline-step-by-step)
  * [Trim reads with trimmomatic](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#1--trim-reads-with-trimmomatic)
  * [Trim the reads further with cutadapt](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#2--trim-the-reads-further-with-cutadapt)
  * [Align reads with bowtie2](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#1--trim-reads-with-trimmomatic)
  * [Filter Alignments](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#1--trim-reads-with-trimmomatic)
  * [Optional step:  Separate subject reads from spike-in reads](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#1--trim-reads-with-trimmomatic)
  * [Make coverage files for genome viewer with deeptools](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#1--trim-reads-with-trimmomatic)
  * [Call peaks with](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#5--call-peaks):
    * [macs2 narrow](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#macs2-narrow)
    * [macs2 broad](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#macs2-broad)
    * [sicer](https://github.com/learn-bioinformatics/ChIPSeq-Cut-and-Run/edit/main/README.md#sicer)

* [Step-by-step instructions on running Snakemake pipeline:](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#step-by-step-instructions-on-running-snakemake-pipeline)
  * [1.  Load slurm and miniconda](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#1--load-slurm-and-miniconda)
  * [2.  Clone repository](https://github.com/SansamLab/Process_HiC_SnakeMake#2--clone-repository)
  * [3.  Start the conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3--start-the-conda-environment)
    * [3A.  FIRST TIME ONLY:  Setup conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3a--first-time-only--setup-conda-environment)
    * [3B.  Activate conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3b--activate-conda-environment)
  * [4.  Modify the job-specific configuration files.](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--modify-the-job-specific-coniguration-files)
    * [4A.  Modify the config/config.yml file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4a--modify-the-configconfigyml-file)
    * [4B.  Modify the config/samples.csv file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4b--modify-the-configsamplescsv-file)
    * [4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4c--if-slurm-resource-changes-are-needed-modify-the-configcluster_configyml-file)
  * [5.  Do a dry run](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--do-a-dry-run)
  * [6.  Make a DAG diagram](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#5--make-a-dag-diagram)
  * [7.  Run on cluster with slurm](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#6--run-on-cluster-with-slurm)

## Requirements

* Packages versions used for tests:
  * trimmomatic/0.35
  * cutadapt/3.2
  * bowtie2/2.3.1
  * samtools/1.14
  * sambamba/0.4.7
  * ucsc/202002225
  * deeptools/3.4.3
* A genome index for bowtie2 [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome). If a spike-in was used, you could index a "hybrid/chimera" genome of the experimental and spike-in species. For example, when running a human cut&run sample with a drosophila genomic DNA spike-in, you could align your reads to an hg19 human genome/dm drosophila genome chimera.
* A path to a .fasta file containing all the adapters, PCR sequences etc. See the [Trimmomatic manual](http://www.usadellab.org/cms/?page=trimmomatic) for more details. [Example For TruSeq Adapters](https://github.com/usadellab/Trimmomatic/tree/main/adapters).
* A path to a .bed file with coordinates of blacklisted and or greylisted regions.
* Pairs of paired-end fastq files for each sample/control. Each ChiP or Cut&Run sample must have a matched control to account for biases in the sequencing process and alignment, such as PCR amplification artifacts, GC biases and alignment artifacts. 

## Run the pipeline step-by-step

### 1.  Trim reads with trimmomatic.

Trimmomatic is performed 1st because it handles quality problems within a sequence, using the sliding window approach. Here we set the sliding window to evaluate 4 base pairs, and check if the average quality within the window is greater than or equal to 15.

```bash
SAMPLE="YOURSAMPLENAMEHERE" `#for the following paired fastq files:  \
                           YOURSAMPLENAMEHERE_R1_001.fastq.gz and 
                           YOURSAMPLENAMEHERE_R2_001.fastq.gz`

trimmomatic PE \
    -threads $THREADS \
    ${SAMPLE}_R1_001.fastq.gz              `# input forward reads`\
    ${SAMPLE}_R2_001.fastq.gz              `# input reverse reads`\
    ${SAMPLE}_R1_trimmed.fastq.gz          `# output trimmed forward reads`\
    ${SAMPLE}_R1_trimmed_orphaned.fastq.gz `# output forward reads that lack corresponding reverse read`\
    ${SAMPLE}_R2_trimmed.fastq.gz          `# output trimmed reverse reads`\
    ${SAMPLE}_R2_trimmed_orphaned.fastq.gz `# ouptut reverse reads that lack corresponding forward reads`\
    ILLUMINACLIP:$ADAPTER_FILE:2:15:4:4:true `# Trim using designated adapter file`\
                                             `# For "2:15:4:4:true":`\
                                             `#      2 is allowed mismatches`\
                                             `#        15 is minimum alignment score for palindromic sequences`\
                                             `#           4 is minimum alignment score between a read and adapter`\
                                             `#             4 is the minimum length of a palindromic "adapter" to be trimmed`\
                                             `#               true means to keep a palindromic read so that the paired files don't get out of sync`\
    LEADING:20          `# Trim 5' bases that have a quality score less than 20`\
    TRAILING:20         `# Trim 3' bases that have a quality score less than 20`\
    SLIDINGWINDOW:4:15  `# Trim anywhere that the average quality score for 4 consecutive bases is less than 15`\
    MINLEN:25           `# Discard reads if shorter than 25 bases`\
```

### 2.  Trim reads further with cutadapt.

The problem is that trimmomatic won't trim adapters matches that are less than 6 base pairs. So, for example, if the adapter is AGATTTACCTAGGAAT and the 3' end of the read is: AATCATATTACAGAT Trimmomatic will NOT trim the AGAT. Cutadapt will allow us to trim these. Note:  An adapter file should be provided. These are likely

```bash
cutadapt \
    --cores=0                                          `# use all available cores`\
    -a file:$ADAPTER_FILE                              `# forward adapter file `\
    -A file:$ADAPTER_FILE                              `# reverse adapter file `\
    -a G{100}                                          `# Helps remove Poly G problem from R1 reads`\
    -A G{100}                                          `# Helps remove Poly G problem from R2 reads`\
    --minimum-length 25                                `# minimum read length after trimming`\
    --quality-cutoff 20,20                             `# threshold for trimming both 5' and 3' `\
    -e 0.2                                             `# mismatch rate 20% (default 0.1) `\
    --output ${SAMPLE}_R1_trimmed2.fastq.gz        `# output forward reads`\
    --paired-output ${SAMPLE}_R2_trimmed2.fastq.gz `# output reverse reads`\
    ${SAMPLE}_R1_trimmed.fastq.gz                      `# input forward reads`\
    ${SAMPLE}_R2_trimmed.fastq.gz                      `# input reverse reads`\
                                                        #*End of cutadapt command*#
```

### 3.  Align reads with bowtie2.
bowtie2 is our aligner. We will map our reads to the provided indexed genome. We then pipe the result from bowtie2 into samtools, with the -b and -S flags. -b will output the result as a BINARY sam file (.bam). Then we run samtools sort, which sorts the reads based on their genomic position.
```bash
bowtie2 \
    --threads "${THREADS}"            `# Number of processing threads to use`\
    --dovetail                        `# Allow overhanging overlaps between R1 and R2 reads (may have resulted from quality trimming)`\
    --phred33                         `# Use the typical quality score scheme`\
    --maxins 2000                     `# Maximum fragment size, including the forward and reverse read and any sequence between them`\
    -x ${BOWTIE2_REF_INDEX}           `# bowtie2 genome index`\
    -1 ${SAMPLE}_R1_trimmed2.fastq.gz `# input forward reads`\
    -2 ${SAMPLE}_R2_trimmed2.fastq.gz `# input reverse reads`\
    |                                 `# pipe`\
samtools view -b -                    `# convert standard input to Binary Alignment Format`\
    |                                 `# pipe`\
samtools sort --threads ${THREADS} - -o ${SAMPLE}.sorted.bam `# Sort BAM file`\
                                       #*End of bowtie2 through samtools piped commands*#
# Index BAM file
samtools index ${SAMPLE}.sorted.bam
```
### 4.  Filter Alignments
Mapping Quality Filter, keep the reads that map with quality >= 10 and removes reads where the paired mate is unmapped, and removes duplicated reads. Output is still .bam format, SAMPLE_NAME.sorted.filt.bam

```bash
sambamba view \
            -h `# include header`\
            -t "${THREADS}" \
            -f bam `# output in BAM format`\
            -F "mapping_quality >= 10 and not mate_is_unmapped and not duplicate" \
            "${SAMPLE}.sorted.bam" > "${SAMPLE}.sorted.filt.bam"

# Index BAM file
samtools index ${SAMPLE}.sorted.filt.bam
```
### Optional step:  Get reads from specific species.
Using awk, we generate the hg19.bed file. This is simply the Start Position (0) to the End Position of each chromosome, in nucleotide bases.
We use it to extract only the human reads from the .bam file, which contains the chimera genome and all the non-human reads.
```
# make human chromosomes bed file from hg19
awk '{print $1 "\t0\t" $2}' /Volumes/shared-refs/hg19/hg19.fasta.fai > temp_hg19_bed
samtools view -b -L temp_hg19_bed ${SAMPLE}.sorted.filt.bam > ${SAMPLE}.sorted.filt.hs.bam

# Index BAM file
samtools index ${SAMPLE}.sorted.filt.hs.bam

rm temp_hg19_bed
```

### 5.  Make coverage files for genome viewer with deeptools.

```bash
BIN_SIZE=10

bamCoverage \
    -b ${SAMPLE}.sorted.filt.hs.bam \
    -o "${SAMPLE}_cpm.bw" \
    --normalizeUsing CPM \
    --binSize $BIN_SIZE \
    --numberOfProcessors max \
    --verbose \
    --blackListFileName $BLACKLIST_FILE \
    --centerReads
```
### 5.  Call peaks.

#### macs2 narrow

```bash
```

#### macs2 broad

```bash
```

#### sicer

```bash
```

## Step-by-step instructions on running Snakemake pipeline:

### 1.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC using slurm with these modules installed.

```bash
ml slurm
ml miniconda
```
### 2.  Clone repository
```bash
git clone https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake.git
# rename folder with project name
mv Cut_And_Run_Analysis_SnakeMake/ My_CutAndRun_Project_Folder/
# change directory into root of your project folder
cd My_CutAndRun_Project_Folder
```
### 3.  Start the conda environment
### 3A.  FIRST TIME ONLY:  Setup conda environment
```bash
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
conda env create -f workflow/envs/CutAndRun_Conda_Environment.yml -p /s/sansam-lab/CutAndRun_Conda_Environment 
```

### 3B.  Activate conda environment
```bash
conda activate /s/sansam-lab/CutAndRun_Conda_Environment
```

### 4.  Modify the job-specific coniguration files.
#### 4A.  Modify the config/config.yml file

You must enter paths to the following:
* bwa_genome:
  * location of bwa indexed genome for the alignment
* chrom_sizes
  * chromosome sizes file
* juicer_RE_file
  * restriction enzyme file generated with juicer

#### 4B.  Modify the config/samples.csv file

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files.

#### 4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file

CPU and memory requests for each rule in the pipeline are detailed in this file. If you are using SLURM, you may need to alter this file to fit your needs/system.

### 5.  Do a dry run.
A dry run produces a text output showing exactly what commands will be executed. Look this over carefully before submitting the full job. It is normal to see warnings about changes made to the code, input, and params.
```bash
snakemake -npr
```

### 6.  Make a DAG diagram.
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### 7.  Run on cluster with slurm.
This snakemake pipeline could be executed without slurm, but if an hpc with slurm is used, the following will start the pipeline with the parameters defined in the config/cluster_config.yml file.
```bash
sbatch --wrap="\
snakemake \
-R \
-j 999 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```

### 8.  Check results, and when finished, exit environment.
The results will be saved to the "results" folder. Look over log files generated in either the logs/ or logs/snakelogs folders (depending on whether slurm was used).
```bash
conda deactivate
```

## Citations

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.


## Credits

## License






This README.md is a work in progress!

TODO:

Get config from a tab-delimited file, which can be read using parallel --colsep='\t' (or comma-delimited).

QUESTIONS:

Do we need a blacklist for spike-ins?
Should we include .bed files for chromosome sizes for human, mouse, and zebrafish UCSC genomes?


NEEDED:

Provenance for blacklist file(s)

Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z
    https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz

Provenance for human-drosophila merged genome files which we (Nic and/or Christopher) created it (i.e. where is(are) our build script(s)?).
