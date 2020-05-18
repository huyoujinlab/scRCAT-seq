#!/bin/bash

fa_file=$1
star_dir=$2
fq_dir=$3
PolyA_db=$4
threshold=$5


## Check whether variables are enough
if [[ ! -n "$fa_file" || ! -n "$star_dir" || ! -n "$fq_dir" || ! -n "$PolyA_db" || ! -n "$threshold" ]]; then    
	echo "Error: no enough variables!
This script is used to call peak, generate features and estimate whether peak is TRUE or FALSE
Usage:
    sh $0 <genome fa file> <STAR index> <fastq file dir> <PolyA_db file> threshold 
    # for example:
    sh callpeak_correction_TSS.sh ~/index/hg38/hg38.fa ~/index/hg38_STAR2.7.3a_index/ ./fastq human.PAS100_hg38.bed 3
    # Caution: file in <fastq file dir> must be fastq type, not fastq.gz! Pair-end files are needed!
    "
    exit
fi





# Data processing for 3' data

#The workflows of data of scCAT-seq 3' and BAT-seq are similar. Here is the scCAT-seq 5' data processing workflow. To see detail imformation of BAT-seq data processing, please see `BAT-seq_3_data_processing.sh`.

## 0. Preparation

#Before process the data, we bulid some directory and move the script to `script_and_log` directory:


#### Create directory
mkdir outdir/
mkdir outdir/three_prime/
mkdir outdir/three_prime/3tail_read_with_tag/
mkdir outdir/three_prime/3tail_read_with_tag_other_strand/
mkdir outdir/three_prime/3tail_read_with_tag_other_strand_withA10_remain_A5
mkdir outdir/three_prime/mapping_outdir/
mkdir outdir/three_prime/final_out/
mkdir outdir/three_prime/annote/
mkdir outdir/three_prime/collapse/


## 1. Find reads with oligo(dT) primer

#Reads with oligo(dT) primer sequence at 5'. We define reads with oligo(dT) primer sequence at 5' as R1 reads. Oligo(dT) primers in scCAT-seq data are listed in `sample_list_tag.txt`:

for i in `ls ${fq_dir}|grep "8N"`
do
        cat ${fq_dir}/${i} | paste - - - - | grep -E $'\t'"GTGGTATCAACGCAGAGT[A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T][A|G|C|T]CTAAGCCTTTTT" | awk -v FS="\t" -v OFS="\n" '{print $1, $2, $3, $4}' > outdir/three_prime/3tail_read_with_tag/${i}_with_tag
done



#Output files are stored in `outdir/three_prime/3tail_read_with_tag/`.

## 2. Find R2 reads

#Perl script `cmpfastq_pe.pl` is used to find R2 reads which its corresponding R1 reads with oligo(dT) primes:


#### Compare R1_with_tag to R2
for i in `ls outdir/three_prime/3tail_read_with_tag/ | grep "R1.fastq_with_tag$"`
do
        perl script/cmpfastq_pe.pl outdir/three_prime/3tail_read_with_tag/${i} ${fq_dir}/${i%R1*}R2.fastq
done

#### Compare R2_with_tag to R1
for i in `ls outdir/three_prime/3tail_read_with_tag/ | grep "R2.fastq_with_tag$"`
do
        perl script/cmpfastq_pe.pl outdir/three_prime/3tail_read_with_tag/${i} ${fq_dir}/${i%R2*}R1.fastq
done

#### Remove useless files
rm outdir/three_prime/3tail_read_with_tag/*out
rm ${fq_dir}/*unique.out
mv ${fq_dir}/*out outdir/three_prime/3tail_read_with_tag_other_strand/


#Output files are stored in `outdir/three_prime/3tail_read_with_tag_other_strand/`.

## 3. Trim A10 at R2 reads

#To trim polyA at 3', we run:


for i in `ls outdir/three_prime/3tail_read_with_tag_other_strand`
do
        python script/find_A10_and_trim.py  -i outdir/three_prime/3tail_read_with_tag_other_strand/${i} -o outdir/three_prime/3tail_read_with_tag_other_strand_withA10_remain_A5/${i}_withA10_remain_A5
done


#Output files are stored in `outdir/three_prime/3tail_read_with_tag_other_strand_withA10_remain_A5/`.

## 4. Mapping

#For Mapping, we run:


for i in `ls outdir/three_prime/3tail_read_with_tag_other_strand_withA10_remain_A5/`
do
        STAR --runThreadN 24 --genomeDir ${star_dir} --genomeLoad LoadAndKeep --readFilesIn outdir/three_prime/3tail_read_with_tag_other_strand_withA10_remain_A5/${i} --outFileNamePrefix outdir/three_prime/mapping_outdir/${i}_ --outSAMtype SAM --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6
done


#Output files are stored in `outdir/three_prime/mapping_outdir/`.


## 5. Remove useless file


for i in `ls outdir/three_prime/mapping_outdir |grep "sam"|grep "R1"`
do
### The first read count
        a=$(wc -l outdir/three_prime/mapping_outdir/${i}|awk '{print $1}')
        echo ${a}

### The second read count
        b=$(wc -l outdir/three_prime/mapping_outdir/${i%%R1*}R2.fastq-common.out_withA10_remain_A5_Aligned.out.sam|awk '{print $1}')
        echo ${b}

### Remove useless files
        if [ ${a} -gt ${b} ]; then
                rm outdir/three_prime/mapping_outdir/${i%%R1*}R2*
        else
                rm outdir/three_prime/mapping_outdir/${i%%R1*}R1*
        fi
done


## 6. Convert SAM to BED

#As `BED` format file can be used as input for `CAGEr` R package, we generate `SAM` to `BED`:


for i in `ls outdir/three_prime/mapping_outdir | grep "sam$"`
do
### Add header and convert to bam
        samtools view -b -T ${fa_file} outdir/three_prime/mapping_outdir/${i} | samtools view -b > outdir/three_prime/final_out/${i}_add_header.bam

### Sort
        samtools sort outdir/three_prime/final_out/${i}_add_header.bam -o outdir/three_prime/final_out/${i}_add_header_sorted.bam

### Build bam index for visualization
        samtools index outdir/three_prime/final_out/${i}_add_header_sorted.bam

### Convert bam into bed
        bedtools bamtobed -i outdir/three_prime/final_out/${i}_add_header_sorted.bam > outdir/three_prime/final_out/${i}_add_header_sorted.bed
done


#Output files are stored in `outdir/three_prime/final_out/`.



## 7. Remove reads mapped to tRNA and rRNA

#Reads that mapped at tRNA and rRNA position are discarded:

for i in `ls  outdir/three_prime/final_out |grep "bed$"`
do
        bedtools subtract -a outdir/three_prime/final_out/${i} -b reference/gencode_hg38_tRNA_rRNA_gene.bed > outdir/three_prime/final_out/${i%.*}_remove_trRNA.bed
done



## 8. find barcode and UMI

for i in `ls outdir/three_prime/final_out | grep "sorted_remove_trRNA.bed$"`
do

        if [ ${i#*.R} == "1.fastq-common.out_withA10_remain_A5_Aligned.out.sam_add_header_sorted_remove_trRNA.bed" ]; then
                python2 script/annotate_UMI_v1.py -N 18 -n 8 -F ${fq_dir}/${i%%.R*}.R2.fastq -ID outdir/three_prime/final_out/${i} -O outdir/three_prime/annote/${i}.annote
        fi
        if [ ${i#*.R} == "2.fastq-common.out_withA10_remain_A5_Aligned.out.sam_add_header_sorted_remove_trRNA.bed" ]; then
                python2 script/annotate_UMI_v1.py -N 18 -n 8 -F ${fq_dir}/${i%%.R*}.R1.fastq -ID outdir/three_prime/final_out/${i} -O outdir/three_prime/annote/${i}.annote
        fi
done


## 9. collapse

for i in `ls outdir/three_prime/annote | grep "annote$"`
do
        Rscript script/collapse_UMI_TES.R  outdir/three_prime/annote/${i}
done

mv outdir/three_prime/annote/*collapse outdir/three_prime/collapse/


## 10. collapse 6 and change name

for i in `ls outdir/three_prime/collapse | grep "collapse$"`
do
        awk '{FS=" "}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' outdir/three_prime/collapse/${i} > outdir/three_prime/collapse/${i}.6
        mv outdir/three_prime/collapse/${i}.6 outdir/three_prime/collapse/${i%%L3*}3_TKD.bed
done








for i in `ls outdir/three_prime/collapse | grep "TKD.bed$"`
do
        bed_file=outdir/three_prime/collapse/${i}
done

filename=${bed_file##*/}
prefix=${filename%%_TKD*}


mkdir outdir/three_prime/peakfile/




## Call peak 
Rscript script/CAGE_dominant_TES.R ${bed_file} ${threshold}


## Extract peak in gene region
bedtools intersect -s -a outdir/three_prime/peakfile/${prefix}_3tail_dominant_tes.bed -b reference/gencode_hg38_all_gene_genebody_and_downstream2k.bed -wa -wb > outdir/three_prime/peakfile/${prefix}_3tail_dominant_tes_genebody_and_downstream2k.bed


## read distribution feature generation
Rscript script/read_distribution_TES.R outdir/three_prime/peakfile/${prefix}_3tail_dominant_tes_genebody_and_downstream2k.bed outdir/three_prime/peakfile/${prefix}_3tail_dominant_tes.bed

echo ${prefix}
## Compare to polydb
bedtools intersect -s -a outdir/three_prime/peakfile/temp_${prefix}_tes.bed -b ${PolyA_db} -wa -wb > outdir/three_prime/peakfile/temp_${prefix}_tes_in_polydb.bed
Rscript script/polydb.R outdir/three_prime/peakfile/tc_${prefix}_3tail.csv


## Add motif information
python script/find_motif_re_TES.py ${fa_file} outdir/three_prime/peakfile/tc_${prefix}_3tail_new.csv outdir/three_prime/peakfile/tc_${prefix}_3tail_new_new.csv


## Unfold motif
Rscript script/unfold_TES.R outdir/three_prime/peakfile/tc_${prefix}_3tail_new_new.csv


## add internal feature
sh script/find_internal_TES.sh outdir/three_prime/peakfile/tc_${prefix}_3tail_final.csv ${fa_file} reference/gencode_hg38_all_gene.bed 

