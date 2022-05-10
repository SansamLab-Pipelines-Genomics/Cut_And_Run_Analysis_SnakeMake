#!/usr/bin/bash -l

# get hg38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

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
zcat GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz dm_ec.fna.gz > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_dm_ec.fna

# index combined genome with bowtie2
## start the bowtie2 aligner....ie "ml bowtie2"
bowtie2-build GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_dm_ec.fna hg38_ecoli_dm_bowtie2Index
