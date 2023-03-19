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

* [Requirements](https://github.com/SansamLab-Pipelines-Genomics/Cut_And_Run_Analysis_SnakeMake/blob/main/docs/requirements.md)
* [Description of individual steps in pipeline](https://github.com/SansamLab-Pipelines-Genomics/Cut_And_Run_Analysis_SnakeMake/blob/main/docs/run_pipeline_stepByStep.md)
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

Provenance for blacklist file(s)

Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z
    https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz
