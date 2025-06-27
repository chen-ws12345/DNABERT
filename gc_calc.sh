#!/bin/bash


bedtools nuc -fi ./genome_ref/T2TCHM13v2_renamed_genomic.fa -bed AR_peaks_VCap.bed > AR_pos_peaks_gc_content.txt

bedtools nuc -fi ./genome_ref/T2TCHM13v2_renamed_genomic.fa -bed random_non_AR_overlapping_regions.bed > random_non_AR_overlapping_regions_gc.txt
