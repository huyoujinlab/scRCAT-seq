#!/bin/bash


fa_file=$1
bed_dir=$2
FANTOM5=$3
threshold=$4

## Check whether variables are enough
if [[ ! -n "$fa_file" || ! -n "$bed_dir" || ! -n "$FANTOM5" || ! -n "$threshold" ]]; then
    echo "Error: no enough variable!
This script is used to call TSS peaks, and generate features which will then be used to predict authentic TSSs.
Usage.
Usage:
    sh $0 <genome fa file> <bed file dir> <FANTOM5 file> threshold
    # for example:
    sh $0 ~/index/hg38/hg38.fa outdir/five_prime/collapse/ reference/tc_hESC.bed 3
    "
    exit
fi

if [ -d "$bed_dir" ]; then
    for i in `ls ${bed_dir} | grep "TKD.bed$"`
    do
        bed_file=${bed_dir}/${i}
    done
fi
if [ -f "$bed_dir" ]; then
    bed_file=${bed_dir}
fi

filename=${bed_file##*/}
prefix=${filename%%_TKD*}
export prefix

mkdir -p outdir/five_prime/peakfile

#######callpeak and feature generation

## Call TSS peaks  
Rscript script/CAGE_dominant_TSS.R ${bed_file} ${threshold}


## Extract peaks related to annotated genes
bedtools intersect -s -a outdir/five_prime/peakfile/${prefix}_5cap_dominant_tss.bed -b reference/gencode_hg38_all_gene_upstream2k_and_genebody.bed -wa -wb > outdir/five_prime/peakfile/${prefix}_5cap_dominant_tss_upstream2k_and_genebody.bed


## Features related to read distribution
Rscript script/read_distribution_TSS.R outdir/five_prime/peakfile/${prefix}_5cap_dominant_tss_upstream2k_and_genebody.bed outdir/five_prime/peakfile/${prefix}_5cap_dominant_tss.bed


## Label the TSS peaks with FANTOM5 DB.
bedtools intersect -s -a outdir/five_prime/peakfile/temp_${prefix}_tss.bed -b ${FANTOM5} -wa -wb > outdir/five_prime/peakfile/temp_${prefix}_tss_in_FANTOM.bed
Rscript script/FANTOM.R outdir/five_prime/peakfile/tc_${prefix}_5cap.csv


## Features related to motif
python script/find_motif_re_TSS.py ${fa_file} outdir/five_prime/peakfile/tc_${prefix}_5cap_new.csv outdir/five_prime/peakfile/tc_${prefix}_5cap_new_new.csv


## Unfold motif
Rscript script/unfold_TSS.R outdir/five_prime/peakfile/tc_${prefix}_5cap_new_new.csv


## Features related to internal TSSs
sh script/find_internal_TSS.sh outdir/five_prime/peakfile/tc_${prefix}_5cap_final.csv ${fa_file} reference/gencode_hg38_all_gene.bed 


