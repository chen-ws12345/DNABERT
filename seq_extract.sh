#!/bin/bash

#bedtools getfasta -fi ./genome_ref/T2TCHM13v2_renamed_genomic.fa -bed AR_peaks_VCap.bed -fo AR_pos_peaks_sequences.fa

#bedtools getfasta -fi ./genome_ref/T2TCHM13v2_renamed_genomic.fa -bed random_non_AR_overlapping_regions_filtered_gc.bed -fo random_non_AR_overlapping_regions_filtered_gc.fa

bedtools getfasta -fi ./genome_ref/T2TCHM13v2_renamed_genomic.fa -bed TRs_w_padding.bed -fo TRs_w_padding_sequences.fa
