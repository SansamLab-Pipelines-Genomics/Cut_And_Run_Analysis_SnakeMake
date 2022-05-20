#!/bin/bash

# if this script is executed from the snakemake results directory it will generate md5sum hashes for important results files

## make directory for the md5sum hash files
mkdir -p md5_hashes

md5sum results/aligned/*.bam > results/md5_hashes/aligned.md5
md5sum results/aligned/*.bai >> results/md5_hashes/aligned.md5
md5sum results/aligned_speciesOfInterest/*.bam > results/md5_hashes/aligned_speciesOfInterest.md5
md5sum results/aligned_speciesOfInterest/*.bai >> results/md5_hashes/aligned_speciesOfInterest.md5
md5sum results/bigwigs_no_spikein/*.bw > results/md5_hashes/bigwigs_no_spikein.md5
md5sum results/bigwigs_spikein/*.bw > results/md5_hashes/bigwigs_spikein.md5
md5sum results/sicer/* > results/md5_hashes/sicer.md5
md5sum results/macs2_normalPeaks/* > results/md5_hashes/macs2_normalPeaks.md5
md5sum results/macs2_broadPeaks/* > results/md5_hashes/macs2_broadPeaks.md5
