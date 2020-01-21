#!/bin/bash

set -e


bed_file=$1
smartseq2_file=$2
wig_file=$3
PolyA_db=$4

filename=${bed_file##*/}
prefix=${filename%%_TKD*}


## Check whether variables are enough
if [[ ! -n "$bed_file" || ! -n "$smartseq2_file" || ! -n "$wig_file" || ! -n "$PolyA_db" ]]; then
    echo "Error: no enough variables!
This script is used to call peak, generate features and estimate whether peak is TRUE or FALSE
Usage:
    sh $0 <bed file> <Smart-seq2 readcount> <Smart-seq2 wig file> <PolyA_db file>
    # for example:
    sh callpeak_correction_TES.sh O41_72_final_3tail.bed O_merge_standard_smart_seq2.count O_merge_standard_smart_seq2_sorted_chr3.wig mouse.PAS100_mm10.bed
    "
    exit
fi



## Call peak 
Rscript CAGE_dominant_TES.R ${bed_file}


## Extract peak in gene region
bedtools intersect -s -a ${prefix}_3tail_dominant_tes.bed -b gencode_mm10_all_gene_genebody_and_downstream2k.bed -wa -wb > ${prefix}_3tail_dominant_tes_genebody_and_downstream2k.bed


## Calculate slope, correaltion, ect.
Rscript cal_slope_TES.R ${prefix}_3tail_dominant_tes_genebody_and_downstream2k.bed ${smartseq2_file} ${wig_file}


## Compare to FANTOM5
bedtools intersect -s -a temp_tes.bed -b ${PolyA_db} -wa -wb > temp_tes_in_polydb.bed
Rscript polydb.R tc_${prefix}_3tail.csv


## Add motif information
python find_motif_re_TES.py ~/index/mm10/mm10.fa tc_${prefix}_3tail_new.csv tc_${prefix}_3tail_new_new.csv


## Unfold motif
Rscript unfold_TES.R tc_${prefix}_3tail_new_new.csv




