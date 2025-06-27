#!/bin/bash

# randomly select 63476 regions across the genome
bedtools random -l 151 -n 10000 -g ./genome_ref/t2t.chrom.sizes > random_regions.bed
# remove regions overlapping with AR peaks
bedtools intersect -a random_regions.bed -b AR_peaks_VCap.bed -v > random_non_AR_overlapping_regions.bed
