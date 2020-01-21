#!/bin/bash
####################################################################################################################
echo
echo
echo
echo "find A5 and remain A5"
echo "########################################################"
echo "########################################################"
echo
echo
echo


mkdir ~/zjw/20190412upload/BAT-seq-dir/trim_and_remainA5_filter_yes/

for i in `ls ~/zjw/20190412upload/BAT-seq-dir/BAT-seq/|grep "fastq$"`
do
python find_A10_and_remain_A5.py -i ~/zjw/20190412upload/BAT-seq-dir/BAT-seq/${i} -o ~/zjw/20190412upload/BAT-seq-dir/trim_and_remainA5_filter_yes/${i}_withA10_remainA5
done
####################################################################################################################


####################################################################################################################
echo
echo
echo
echo "mapping"
echo "########################################################"
echo "########################################################"
echo
echo
echo


mkdir ~/zjw/20190412upload/BAT-seq-dir/mapping_output_filter_yes/
for i in `ls ~/zjw/20190412upload/BAT-seq-dir/trim_and_remainA5_filter_yes/`
do
STAR --runThreadN 24 --genomeDir ~/index/mm10_ERCC92trimpolyA_STAR --genomeLoad LoadAndKeep --readFilesIn ~/zjw/20190412upload/BAT-seq-dir/trim_and_remainA5_filter_yes/${i} --outFileNamePrefix ~/zjw/20190412upload/BAT-seq-dir/mapping_output_filter_yes/${i}_ --outSAMtype SAM --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6--outFilterMatchNminOverLread 0.6
done
####################################################################################################################


####################################################################################################################
echo
echo
echo
echo "extract uniquely map"
echo "########################################################"
echo "########################################################"
echo
echo
echo


mkdir ~/zjw/20190412upload/BAT-seq-dir/extract_uniquely_map_filter_yes/
for i in `ls ~/zjw/20190412upload/BAT-seq-dir/mapping_output_filter_yes/ | grep "Aligned.out.sam"`
do
samtools view ~/zjw/20190412upload/BAT-seq-dir/mapping_output_filter_yes/${i} | grep 'NH:i:1'$'\t''' > ~/zjw/20190412upload/BAT-seq-dir/extract_uniquely_map_filter_yes/${i}_extract_uniquely_map.sam
done
####################################################################################################################



####################################################################################################################
echo
echo
echo
echo "split plus minus"
echo "########################################################"
echo "########################################################"
echo
echo
echo


mkdir ~/zjw/20190412upload/BAT-seq-dir/split_plus_minus_filter_yes/
for i in `ls ~/zjw/20190412upload/BAT-seq-dir/extract_uniquely_map_filter_yes/ | grep "sam"`
do
cat ~/zjw/20190412upload/BAT-seq-dir/extract_uniquely_map_filter_yes/${i} | awk '{FS=" "}{if ($2==0 || $2==256){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' > ~/zjw/20190412upload/BAT-seq-dir/split_plus_minus_filter_yes/${i}_plus

cat ~/zjw/20190412upload/BAT-seq-dir/extract_uniquely_map_filter_yes/${i} | awk '{FS=" "}{if ($2==16 || $2==272){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' > ~/zjw/20190412upload/BAT-seq-dir/split_plus_minus_filter_yes/${i}_minus
done
####################################################################################################################



####################################################################################################################
echo
echo
echo
echo "extract mismatch"
echo "########################################################"
echo "########################################################"
echo
echo
echo


mkdir ~/zjw/20190412upload/BAT-seq-dir/extract_mismatch_filter_yes/

for i in `ls ~/zjw/20190412upload/BAT-seq-dir/split_plus_minus_filter_yes/ | grep "extract_uniquely_map.sam_plus"`
do
python extractmismatch_plus.py -i ~/zjw/20190412upload/BAT-seq-dir/split_plus_minus_filter_yes/${i} -o ~/zjw/20190412upload/BAT-seq-dir/extract_mismatch_filter_yes/${i}_extractmismatch

python extractmismatch_minus.py -i ~/zjw/20190412upload/BAT-seq-dir/split_plus_minus_filter_yes/${i%_*}_minus -o ~/zjw/20190412upload/BAT-seq-dir/extract_mismatch_filter_yes/${i%_*}_minus_extractmismatch

cat ~/zjw/20190412upload/BAT-seq-dir/extract_mismatch_filter_yes/${i}_extractmismatch ~/zjw/20190412upload/BAT-seq-dir/extract_mismatch_filter_yes/${i%_*}_minus_extractmismatch > ~/zjw/20190412upload/BAT-seq-dir/extract_mismatch_filter_yes/${i%_*}_extractmismatch
done
####################################################################################################################


####################################################################################################################
echo
echo
echo
echo "final out"
echo "########################################################"
echo "########################################################"
echo
echo
echo


mkdir ~/zjw/20190412upload/BAT-seq-dir/final_out_filter_yes/
for i in `ls ~/zjw/20190412upload/BAT-seq-dir/extract_mismatch_filter_yes/ | grep "sam_extractmismatch"`
do
samtools view -b -T ~/index/mm10_ERCC92/mm10_ERCC92.fa ~/zjw/20190412upload/BAT-seq-dir/extract_mismatch_filter_yes/${i} | samtools view -b > ~/zjw/20190412upload/BAT-seq-dir/final_out_filter_yes/${i}_add_header.bam

samtools sort ~/zjw/20190412upload/BAT-seq-dir/final_out_filter_yes/${i}_add_header.bam -o ~/zjw/20190412upload/BAT-seq-dir/final_out_filter_yes/${i}_add_header_sorted.bam

samtools index ~/zjw/20190412upload/BAT-seq-dir/final_out_filter_yes/${i}_add_header_sorted.bam

bedtools bamtobed -i ~/zjw/20190412upload/BAT-seq-dir/final_out_filter_yes/${i}_add_header_sorted.bam > ~/zjw/20190412upload/BAT-seq-dir/final_out_filter_yes/${i}_add_header_sorted.bed
done

####################################################################################################################
