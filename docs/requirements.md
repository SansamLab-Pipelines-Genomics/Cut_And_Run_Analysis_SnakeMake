
## Requirements

* Packages versions used for tests:
  * trimmomatic/0.35
  * cutadapt/3.2
  * bowtie2/2.3.1
  * samtools/1.14
  * sambamba/0.4.7
  * ucsc/202002225
  * deeptools/3.4.3
  * fastqc/0.11.9
* A genome index for bowtie2 [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome). If a spike-in was used, you could index a "hybrid/chimera" genome of the experimental and spike-in species. For example, when running a human cut&run sample with a drosophila genomic DNA spike-in, you could align your reads to an hg19 human genome/dm drosophila genome chimera.
* A path to a .fasta file containing all the adapters, PCR sequences etc. See the [Trimmomatic manual](http://www.usadellab.org/cms/?page=trimmomatic) for more details. [Example For TruSeq Adapters](https://github.com/usadellab/Trimmomatic/tree/main/adapters).
* A path to a .bed file with coordinates of blacklisted and or greylisted regions.
* Pairs of paired-end fastq files for each sample/control. Each ChiP or Cut&Run sample must have a matched control to account for biases in the sequencing process and alignment, such as PCR amplification artifacts, GC biases and alignment artifacts. 
