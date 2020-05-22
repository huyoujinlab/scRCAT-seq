#!/bin/bash



fa_file=$1
bed_dir=$2
PolyA_db=$3
threshold=$4

## Check variables
if [[ ! -n "$fa_file" || ! -n "$bed_dir" || ! -n "$PolyA_db" || ! -n "$threshold" ]]; then    
	echo "Error: no enough variables!
This script is used to call peak, and generate features which will be used to  build machine learning models to predict authentic TES peaks.
Usage:
    sh $0 <genome fa file> <bed file dir> <PolyA_db file> threshold 
    # for example:
    sh callpeak_correction_TSS.sh ~/index/hg38/hg38.fa outdir/three_prime/collapse/ reference/human.PAS100_hg38.bed 3
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

mkdir -p outdir/three_prime/peakfile/




## Call TES peaks 
Rscript script/CAGE_dominant_TES.R ${bed_file} ${threshold}


## Extract TES peaks located related to annotated genes
bedtools intersect -s -a outdir/three_prime/peakfile/${prefix}_3tail_dominant_tes.bed -b reference/gencode_hg38_all_gene_genebody_and_downstream2k.bed -wa -wb > outdir/three_prime/peakfile/${prefix}_3tail_dominant_tes_genebody_and_downstream2k.bed



##Generate features related to read distribution 
Rscript script/read_distribution_TES.R outdir/three_prime/peakfile/${prefix}_3tail_dominant_tes_genebody_and_downstream2k.bed outdir/three_prime/peakfile/${prefix}_3tail_dominant_tes.bed

echo ${prefix}
## Label the peaks with polydb
bedtools intersect -s -a outdir/three_prime/peakfile/temp_${prefix}_tes.bed -b ${PolyA_db} -wa -wb > outdir/three_prime/peakfile/temp_${prefix}_tes_in_polydb.bed
Rscript script/polydb.R outdir/three_prime/peakfile/tc_${prefix}_3tail.csv


## Generate features related to motif around the peaks
python script/find_motif_re_TES.py ${fa_file} outdir/three_prime/peakfile/tc_${prefix}_3tail_new.csv outdir/three_prime/peakfile/tc_${prefix}_3tail_new_new.csv


## Unfold motif
Rscript script/unfold_TES.R outdir/three_prime/peakfile/tc_${prefix}_3tail_new_new.csv


## generate features related to internal priming events
sh script/find_internal_TES.sh outdir/three_prime/peakfile/tc_${prefix}_3tail_final.csv ${fa_file} reference/gencode_hg38_all_gene.bed 

