# Build and Index a Multi-genome

## Purpose
To concatenate your genome-of-interest with spike-in genomes for Cut&Run alignments.

## Scripts
### build_hg38_ec_dm_multi_genome.sh
#### Genomes:
* Human:  GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set
* Drosophila:  GCF_000001215.4_Release_6_plus_ISO1_MT_genomic
* E. coli:  GCF_000008865.2_ASM886v2_genomic
#### Prefixes added to chromosomes:
This following prefixes are added to the chromosomes to parse reads aligning to each genome.
* Human:  None
* Drosophila:  "DM_"
* E. coli:  "EC_"

#### Executing on an hpc with slurm
```bash
ml bowtie2
sbatch --cpus-per-task 8 --mem 64Gb build_hg38_ec_dm_multi_genome.sh
```
### build_hg38_ec_multi_genome.sh
#### Genomes:
* Human:  GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set
* E. coli:  GCF_000008865.2_ASM886v2_genomic
#### Prefixes added to chromosomes:
This following prefixes are added to the chromosomes to parse reads aligning to each genome.
* Human:  None
* E. coli:  "EC_"

#### Executing on an hpc with slurm
```bash
ml bowtie2
sbatch --cpus-per-task 8 --mem 64Gb build_hg38_ec_multi_genome.sh
```
