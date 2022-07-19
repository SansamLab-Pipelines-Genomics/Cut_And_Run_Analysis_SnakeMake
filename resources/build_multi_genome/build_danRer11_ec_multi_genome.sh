#!/usr/bin/bash -l

# get hg38
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Danio_rerio/latest_assembly_versions/GCA_000002035.4_GRCz11/GCA_000002035.4_GRCz11_genomic.fna.gz

# get e coli genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
# modify chromosome names
zcat GCF_000008865.2_ASM886v2_genomic.fna.gz | sed 's/>/>EC_/g' > ec.fasta

# concatenate fasta files
gzip > ec.fna.gz
zcat GCA_000002035.4_GRCz11_genomic.fna.gz ec.fna.gz > GCA_000002035.4_GRCz11_genomic.fna

# make directory for index
mkdir danRerio_ecoli_bowtie2Index
cd danRerio_ecoli_bowtie2Index

# index combined genome with bowtie2
## start the bowtie2 aligner....ie "ml bowtie2"
bowtie2-build --threads 8 ../GCA_000002035.4_GRCz11_genomic.fna danRerio_ecoli_bowtie2Index
