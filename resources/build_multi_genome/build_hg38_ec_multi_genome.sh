#!/usr/bin/bash -l

# get hg38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

# get e coli genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
# modify chromosome names
zcat GCF_000008865.2_ASM886v2_genomic.fna.gz | sed 's/>/>EC_/g' | gzip > ec.fna.gz

# concatenate fasta files
zcat GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz ec.fna.gz > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_ec.fna

# make directory for index
mkdir hg38_ecoli_bowtie2Index
cd hg38_ecoli_bowtie2Index

# index combined genome with bowtie2
## start the bowtie2 aligner....ie "ml bowtie2"
bowtie2-build --threads 8 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_ec.fna hg38_ecoli_bowtie2Index
