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
