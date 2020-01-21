#!/bin/bash

set -e


bed_file=$1
smartseq2_file=$2
wig_file=$3
FANTOM5=$4

filename=${bed_file##*/}
prefix=${filename%%_TKD*}


## Check whether variables are enough
if [[ ! -n "$bed_file" || ! -n "$smartseq2_file" || ! -n "$wig_file" || ! -n "$FANTOM5" ]]; then
    echo "Error: no enough variables!
This script is used to call peak, generate features and estimate whether peak is TRUE or FALSE
Usage:
    sh $0 <bed file> <Smart-seq2 readcount> <Smart-seq2 wig file> <FANTOM5 file>
    # for example:
    sh callpeak_correction_TSS.sh O41_72_final_5cap.bed O_merge_standard_smart_seq2.count O_merge_standard_smart_seq2_sorted_chr3.wig tc_ovary.bed
    "
    exit
fi



## Call peak  
Rscript CAGE_dominant_TSS.R ${bed_file}


## Extract peak in gene region
bedtools intersect -s -a ${prefix}_5cap_dominant_tss.bed -b gencode_mm10_all_gene_ustream2k_and_genebody.bed -wa -wb > ${prefix}_5cap_dominant_tss_upstream2k_and_genebody.bed


## Calculate slope, correaltion, ect.
Rscript cal_slope_TSS.R ${prefix}_5cap_dominant_tss_upstream2k_and_genebody.bed ${smartseq2_file} ${wig_file}


## Compare to FANTOM5
bedtools intersect -s -a temp_tss.bed -b ${FANTOM5} -wa -wb > temp_tss_in_FANTOM.bed
Rscript FANTOM.R tc_${prefix}_5cap.csv


## Add motif information
python find_motif_re_TSS.py ~/index/mm10/mm10.fa tc_${prefix}_5cap_new.csv tc_${prefix}_5cap_new_new.csv


## Unfold motif
Rscript unfold_TSS.R tc_${prefix}_5cap_new_new.csv




