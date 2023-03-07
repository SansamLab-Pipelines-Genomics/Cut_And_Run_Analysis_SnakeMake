![Release](https://img.shields.io/github/v/release/SansamLab/Cut_And_Run_Analysis_SnakeMake?include_prereleases)
![ReleaseDate](https://img.shields.io/github/release-date/SansamLab/Cut_And_Run_Analysis_SnakeMake)
![Size](https://img.shields.io/github/repo-size/SansamLab/Cut_And_Run_Analysis_SnakeMake)
![License](https://img.shields.io/github/license/SansamLab/Cut_And_Run_Analysis_SnakeMake)
![LastCommit](https://img.shields.io/github/last-commit/SansamLab/Cut_And_Run_Analysis_SnakeMake)
![Downloads](https://img.shields.io/github/downloads/SansamLab/Cut_And_Run_Analysis_SnakeMake/total)
![OpenIssues](https://img.shields.io/github/issues-raw/SansamLab/Cut_And_Run_Analysis_SnakeMake)
[![DOI](https://zenodo.org/badge/468099411.svg)](https://zenodo.org/badge/latestdoi/468099411)


# Cut&Run Pipeline

<figure>
<img src="graphic.png" alt="Trulli" style="width:60%">
<figcaption align = "bottom" style=" text-align : right">image made in BioRender.com</figcaption>
</figure>

## Project Description
Cut_And_Run_Analysis_SnakeMake processes short whole-genome sequencing reads from Cut&Run. The pipeline generates trimmed fastq files, genome alignments, coverage files, and peak calls. To enable step-by-step data processing, we describe each of the individual data processing steps. We provide a Snakemake pipeline with clearly defined dependencies and Anaconda environments to automate the pipeline. We include a compact dataset in the repository that you may use to test the pipeline. We also provide an example detailing how Cut_And_Run_Analysis_SnakeMake can be used to process publicly available Cut&Run data. The original workflow was provided courtesy of the Sansam Lab.

## Table of Contents

* [Requirements](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#requirements)
* [Description of individual steps in pipeline](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#run-the-pipeline-step-by-step)
  * [Trim reads with trimmomatic](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#1--trim-reads-with-trimmomatic)
  * [Trim the reads further with cutadapt](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#2--trim-reads-further-with-cutadapt)
  * [Align reads with bowtie2](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#3--align-reads-with-bowtie2)
  * [Filter Alignments](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#4--filter-alignments)
  * [Optional step:  Separate subject reads from spike-in reads](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#optional-step--get-reads-from-specific-species)
  * [Make coverage files for genome viewer with deeptools](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#5--make-coverage-files-for-genome-viewer-with-deeptools)
  * [Call peaks with](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#5--call-peaks):
    * [macs2 narrow](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#macs2-narrow)
    * [macs2 broad](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#macs2-broad)
    * [sicer](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#sicer)

* [Step-by-step instructions on running Snakemake pipeline:](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#step-by-step-instructions-on-running-snakemake-pipeline)
  * [1.  Load slurm and miniconda](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#1--load-slurm-and-miniconda)
  * [2.  Clone repository](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#2--clone-repository)
  * [3.  Start the conda environment](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#3--start-the-conda-environment)
    * [3A.  FIRST TIME ONLY:  Setup conda environment](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#3a--first-time-only--setup-conda-environment)
    * [3B.  Activate conda environment](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#3b--activate-conda-environment)
  * [4.  Modify the job-specific configuration files.](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#4--modify-the-job-specific-configuration-files)
    * [4A.  Modify the config/config.yml file](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#4a--modify-the-configconfigyml-file)
    * [4B.  Modify the config/samples.csv file](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#4b--modify-the-configsamplescsv-file)
    * [4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#4c--if-slurm-resource-changes-are-needed-modify-the-configcluster_configyml-file)
  * [5.  Do a dry run](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#5--do-a-dry-run)
  * [6.  Make a DAG diagram](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#6--make-a-dag-diagram)
  * [7.  Run on cluster with slurm](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#7--run-on-cluster-with-slurm)
  * [8.  Check results, and when finished, exit environment](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#8--check-results-and-when-finished-exit-environment)
* [Citations](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#citations)
* [License](https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake#license)

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
    -basein ${SAMPLE}_R1_001.fastq.gz        `# input read 1 filename (because of -basein, it will also use read 2) `\
    -baseout ${SAMPLE}_trimmed.fastq.gz      `# base output name for trimmed reads (1P = paired read 1, 1U = unpaired read 1, etc)`\
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
macs2 callpeak \
    -t {input.treatment} \
    -c {input.control} \
    -f BAMPE \
    -g {params.effective_genome_size} \
    -n {params.sample_name}"_"{params.minimum_FDR_cutoff} \
    -q {params.minimum_FDR_cutoff} \
    --outdir results/macs2_normalPeaks/
```

#### macs2 broad

```bash
macs2 callpeak \
    -t {input.treatment} \
    -c {input.control} \
    -f BAMPE \
    -g {params.effective_genome_size} \
    -n {params.sample_name}"_"{params.minimum_FDR_cutoff} \
    -q {params.minimum_FDR_cutoff} \
    --broad \
    --outdir results/macs2_broadPeaks/
```

#### sicer

```bash
sicer \
  -t {input.treatment} \
  -c {input.control} \
  -s {params.sicer_genome} \
  -w {params.sicer_windowSize} \
  -f {params.sicer_fragmentSize} \
  -fdr {params.sicer_fdr} \
  -o results/sicer/ \
  -g {params.sicer_gapSize} \
  -cpu 12
```

# Step-by-step instructions on running Snakemake pipeline with ***Conda Environments***:

### 1.  Load slurm, python, pandas, and numpy
Note. The commands to do this will be different on your machine. These commands are specific to an HPC using slurm with these modules installed.
```bash
# make sure that there are no modules currently loaded
module purge
# load modules to run snakemake file (snakemake was installed in python on this specific cluster)
module load slurm python/3.7.0  pandas/1.0.3  numpy/1.18.2
```

### 2.  Make working directory
```bash
# can name this anything you want
mkdir CutandRunYY1_CondaEnv
# change to the working directory
cd CutandRunYY1_CondaEnv
```

### 3.  Clone repository
```bash
# once you are in the directory you would like to work on the project clone the github repo
git clone https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake.git
# Step below to rename the folder are not necessary but can help organize your files
# rename folder with project name
mv Cut_And_Run_Analysis_SnakeMake/ My_CutAndRun_Project_Folder/
# change directory into root of your project folder
# You will need to be in this directory to run the snakefile
cd My_CutAndRun_Project_Folder
```
** Note. Jump to 5 if conda enviornment is already created **

### 4.  FIRST TIME ONLY:
### 4.  Setup conda environment
Note. This will set up the conda evniornments in what ever location you specify in the --conda-prefix. This step only needs to be done once. After you have this enviornment folder setup you only need to specify its location in step 5.
```bash
# set up conda envs (only run the first time you set up)
sbatch --mem 32G \
--wrap="\
snakemake \
--cores all \
--use-conda \
--conda-prefix /Volumes/Sansam/hpc-nobackup/condaEnvs/CutandRun_KB \
--conda-create-envs-only \
--conda-frontend conda"
```

### 5.  Modify the job-specific configuration files.
#### 5A.  Modify the config/config.yml file

You must enter paths to the following:
* bowtie2_genome:
  * location of bwa indexed genome for the alignment
* blacklist file
  * path to blacklist file to be used for peak calling with macs2 and for making coverage files
* trimmomatic_adapterfile:
  * location of adapter file for trimmomatic (included with trimmomatic for truseq adapters)
* cutadapt_adapterfile:
  * location of adapter file for cutadapt (you may use the same file for trimmomatic)
* sicer:
  * sicer_genome:
  * sicer window size:
  * sicer fragment size:
  * sicer fdr:
  * sicer_gapSize:

#### 5B.  Modify the config/samples.csv file
Note. Make sure to rename sample file by removing "_template"

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files. Columns 2 and 3 are the paths to the paired fastq files. Column 4 is the sample type (either "treatment" or "control"). Column 5 is the name of the corresponding Control sample for each treated sample (use "NA" if the sample is a control). Finally, in the last column put the same sample name you used in the first column to label your sample. This last column is only used when the same sample has more reads in another fastq. In that case you use the sample name to indicate which reads to merge. In the example below, the two treated samples share the same control but, it is possible to include multple controls if that fits your experimental design.

| sample      | fastq1              | fastq2              | sampleType | Control   | merged_sample |
|-------------|---------------------|---------------------|------------|-----------|---------------|
| testSample  | sample_R1.fastq.gz  | sample_R2.fastq.gz  | treatment  | testInput | testSample    |
| testSample2 | sample2_R1.fastq.gz | sample2_R2.fastq.gz | treatment  | testInput | testSample2   |
| testInput   | input_R1.fastq.gz   | input_R2.fastq.gz   | control    | NA        | testInput     |

#### 5C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file

CPU and memory requests for each rule in the pipeline are detailed in this file. If you are using SLURM, you may need to alter this file to fit your needs/system.

### 6.  Do a dry run
A dry run produces a text output showing exactly what commands will be executed. Look this over carefully before submitting the full job. It is normal to see warnings about changes made to the code, input, and params
```bash
snakemake -npr
```

### 7.  Make a DAG diagram
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### 8.  Run on cluster with slurm
This snakemake pipeline could be executed without slurm, but if an hpc with slurm is used, the following will start the pipeline with the parameters defined in the config/cluster_config.yml file.
Note. Make sure to check the --conda-prefix location. This needs to match the location you used when setting up the conda envs
```bash
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-conda \
--conda-prefix /Volumes/Sansam/hpc-nobackup/condaEnvs/CutandRun_KB \
--conda-frontend conda \
--latency-wait 100 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```

### 9.  Check results
The results will be saved to the "results" folder. Look over log files generated in either the logs/ or logs/snakelogs folders (depending on whether slurm was used)
Note. There needs to be an empty results folder in already in your working folder or it will not populate with the files being generated. If this is your first run it will already be there from the git clone step 





# Step-by-step instructions on running Snakemake pipeline with ***Environment Modules***:

### 1.  Make working directory and change to that directory
```bash
# make workingdirectory
mkdir CutandRunYY1_Envmod
# change to working directory
cd CutandRunYY1_Envmod
```

### 2.  Load modules needed for snakemake
```bash
# make sure no modules are already loaded
module purge
# load moduels to run snakemake
module load slurm python/3.7.0  pandas/1.0.3  numpy/1.18.2
```

### 3.  Clone repository
```bash
# once you are in the directory you would like to work on the project clone the github repo
git clone https://github.com/SansamLab/Cut_And_Run_Analysis_SnakeMake.git
# Step below to rename the folder are not necessary but can help organize your files
# rename folder with project name
mv Cut_And_Run_Analysis_SnakeMake/ My_CutAndRun_Project_Folder/
# change directory into root of your project folder
# You will need to be in this directory to run the snakefile
cd My_CutAndRun_Project_Folder
```

### 4.  Modify the job-specific configuration files.
#### 4A.  Modify the config/config.yml file

You must enter paths to the following:
* bowtie2_genome:
  * location of bwa indexed genome for the alignment
* blacklist file
  * path to blacklist file to be used for peak calling with macs2 and for making coverage files
* trimmomatic_adapterfile:
  * location of adapter file for trimmomatic (included with trimmomatic for truseq adapters)
* cutadapt_adapterfile:
  * location of adapter file for cutadapt (you may use the same file for trimmomatic)
* sicer:
  * sicer_genome:
  * sicer window size:
  * sicer fragment size:
  * sicer fdr:
  * sicer_gapSize:

#### 4B.  Modify the config/samples.csv file
Note. Make sure to rename sample file by removing "_template"

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files. Columns 2 and 3 are the paths to the paired fastq files. Column 4 is the sample type (either "treatment" or "control"). Column 5 is the name of the corresponding Control sample for each treated sample (use "NA" if the sample is a control). Finally, in the last column put the same sample name you used in the first column to label your sample. This last column is only used when the same sample has more reads in another fastq. In that case you use the sample name to indicate which reads to merge. In the example below, the two treated samples share the same control but, it is possible to include multple controls if that fits your experimental design.

| sample      | fastq1              | fastq2              | sampleType | Control   | merged_sample |
|-------------|---------------------|---------------------|------------|-----------|---------------|
| testSample  | sample_R1.fastq.gz  | sample_R2.fastq.gz  | treatment  | testInput | testSample    |
| testSample2 | sample2_R1.fastq.gz | sample2_R2.fastq.gz | treatment  | testInput | testSample2   |
| testInput   | input_R1.fastq.gz   | input_R2.fastq.gz   | control    | NA        | testInput     |

#### 4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file

CPU and memory requests for each rule in the pipeline are detailed in this file. If you are using SLURM, you may need to alter this file to fit your needs/system.

### 5.  Do a dry run
A dry run produces a text output showing exactly what commands will be executed. Look this over carefully before submitting the full job. It is normal to see warnings about changes made to the code, input, and params
```bash
snakemake -npr
```

### 6.  Make a DAG diagram
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### 7.  Run on cluster with slurm
```
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-envmodules \
--latency-wait 100 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```

### 8.  Check results
The results will be saved to the "results" folder. Look over log files generated in either the logs/ or logs/snakelogs folders (depending on whether slurm was used)
Note. There needs to be an empty results folder in already in your working folder or it will not populate with the files being generated. If this is your first run it will already be there from the git clone step 



## Citations

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal, 17(1), 10. https://doi.org/10.14806/ej.17.1.200

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of samtools and bcftools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

Tarasov, A., Vilella, A. J., Cuppen, E., Nijman, I. J., & Prins, P. (2015). Sambamba: Fast processing of NGS alignment formats. Bioinformatics, 31(12), 2032–2034. https://doi.org/10.1093/bioinformatics/btv098

Ramírez, F., Ryan, D. P., Grüning, B., Bhardwaj, V., Kilpert, F., Richter, A. S., Heyne, S., Dündar, F., & Manke, T. (2016). Deeptools2: A next generation web server for deep-sequencing data analysis. Nucleic Acids Research, 44(W1), W160–W165. https://doi.org/10.1093/nar/gkw257

Zang, C., Schones, D. E., Zeng, C., Cui, K., Zhao, K., & Peng, W. (2009). A clustering approach for identification of enriched domains from histone modification ChIP-Seq data. Bioinformatics, 25(15), 1952–1958. https://doi.org/10.1093/bioinformatics/btp340

Zhang, Y., Liu, T., Meyer, C. A., Eeckhoute, J., Johnson, D. S., Bernstein, B. E., Nusbaum, C., Myers, R. M., Brown, M., Li, W., & Liu, X. S. (2008). Model-based analysis of chip-seq(Macs). Genome Biology, 9(9), R137. https://doi.org/10.1186/gb-2008-9-9-r137

Amemiya, H. M., Kundaje, A., & Boyle, A. P. (2019). The encode blacklist: Identification of problematic regions of the genome. Scientific Reports, 9(1), 9354. https://doi.org/10.1038/s41598-019-45839-z

Boyle, A. (2018). Boyle-lab/blacklist: Official encode blacklist release for publication. Zenodo. https://doi.org/10.5281/ZENODO.1491733

Köster, J., & Rahmann, S. (2012). Snakemake—A scalable bioinformatics workflow engine. Bioinformatics (Oxford, England), 28(19), 2520–2522. https://doi.org/10.1093/bioinformatics/bts480

Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. Retrieved from https://docs.anaconda.com/

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
