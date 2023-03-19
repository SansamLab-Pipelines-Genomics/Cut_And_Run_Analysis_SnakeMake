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
