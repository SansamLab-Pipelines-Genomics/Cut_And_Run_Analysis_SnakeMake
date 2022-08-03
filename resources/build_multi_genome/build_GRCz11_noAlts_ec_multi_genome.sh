#!/usr/bin/bash -l

#### 1.1 Load necessary software
ml seqtk bwa python


#### 1.2 Copy GRCz11 fasta to local directory
sbatch --wrap=\
"wget -e robots=off --reject 'index.html' https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz"


#### 1.3. Get primary sequences (non-Alts) and rename chromosomes
sbatch --wrap=\
"zcat GCF_000002035.6_GRCz11_genomic.fna.gz | seqtk seq -l 3000000000 | grep Primary -A 1 | sed 's/>/>chrUn/g' | seqtk seq -l 100 | sed 's/^.*chromosome />chr/g' | sed 's/,.*$//g' | gzip > GCF_000002035.6_GRCz11_primary_genomic.fna.gz"


############################################################################


# get e coli genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
# modify chromosome names
zcat GCF_000008865.2_ASM886v2_genomic.fna.gz | sed 's/>/>EC_/g' | gzip > ec.fna.gz

# concatenate fasta files
zcat GCF_000002035.6_GRCz11_primary_genomic.fna.gz ec.fna.gz > GCF_000002035.6_GRCz11_primary_genomic_ec.fna

# make directory for index
mkdir GRCz11_ecoli_bowtie2Index
cd GRCz11_ecoli_bowtie2Index

# index combined genome with bowtie2
## start the bowtie2 aligner....ie "ml bowtie2"
bowtie2-build --threads 8 ../GCF_000002035.6_GRCz11_primary_genomic_ec.fna GRCz11_ecoli_bowtie2Index