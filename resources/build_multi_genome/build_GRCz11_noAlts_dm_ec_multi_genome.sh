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


# get drosophila genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
# modify chromosome names
zcat GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | sed 's/>/>DM_/g' > dm.fasta

# get e coli genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
# modify chromosome names
zcat GCF_000008865.2_ASM886v2_genomic.fna.gz | sed 's/>/>EC_/g' > ec.fasta

# concatenate fasta files
cat dm.fasta ec.fasta | gzip > dm_ec.fna.gz
zcat GCF_000002035.6_GRCz11_primary_genomic.fna.gz dm_ec.fna.gz > GCF_000002035.6_GRCz11_primary_genomic_dm_ec.fna

# make directory for index
mkdir GRCz11_dm_ecoli_bowtie2Index
cd GRCz11_dm_ecoli_bowtie2Index

# index combined genome with bowtie2
## start the bowtie2 aligner....ie "ml bowtie2"
bowtie2-build --threads 8 ../GCF_000002035.6_GRCz11_primary_genomic_dm_ec.fna GRCz11_dm_ecoli_bowtie2Index