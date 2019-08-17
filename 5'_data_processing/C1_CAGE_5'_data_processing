#!/bin/bash

echo "yes filter"
####################################################################################################################
echo
echo
echo
echo "trim UMI"
echo "########################################################"
echo "########################################################"
echo
echo
echo

mkdir ~/zjw/20190416upload/C1_CAGE-dir/trim_UMI_filter_yes/
for i in `ls  ~/zjw/20190416upload/C1_CAGE-dir/C1_CAGE/ | grep "fastq$"`
do
cutadapt -u 8 -o ~/zjw/20190416upload/C1_CAGE-dir/trim_UMI_filter_yes/${i}.trim_UMI  ~/zjw/20190416upload/C1_CAGE-dir/C1_CAGE/${i}
done
####################################################################################################################



####################################################################################################################
echo
echo
echo
echo "find_TATAGGG"
echo "########################################################"
echo "########################################################"
echo
echo
echo
mkdir ~/zjw/20190416upload/C1_CAGE-dir/find_TATAGGG_filter_yes/
for i in `ls  ~/zjw/20190416upload/C1_CAGE-dir/trim_UMI_filter_yes/`
do
cat ~/zjw/20190416upload/C1_CAGE-dir/trim_UMI_filter_yes/${i} | paste - - - - | grep -E $'\t''TATAGGG' | awk -v FS="\t" -v OFS="\n" '{print $1, $2, $3, $4}' > ~/zjw/20190416upload/C1_CAGE-dir/find_TATAGGG_filter_yes/${i}_with_TATAGGG.fastq
done
####################################################################################################################



####################################################################################################################
echo
echo
echo
echo "trim TATA"
echo "########################################################"
echo "########################################################"
echo
echo
echo

mkdir ~/zjw/20190416upload/C1_CAGE-dir/trim_TATA_filter_yes/
for i in `ls  ~/zjw/20190416upload/C1_CAGE-dir/find_TATAGGG_filter_yes/ | grep "fastq$"`
do
cutadapt -u 4 -o ~/zjw/20190416upload/C1_CAGE-dir/trim_TATA_filter_yes/${i}.trim_UMI  ~/zjw/20190416upload/C1_CAGE-dir/find_TATAGGG_filter_yes/${i}
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


mkdir ~/zjw/20190416upload/C1_CAGE-dir/mapping_output_filter_yes/
for i in `ls ~/zjw/20190416upload/C1_CAGE-dir/trim_TATA_filter_yes/ | grep "trim"`
do
STAR --runThreadN 24 --genomeDir ~/index/mm10_ERCC92_STAR --genomeLoad LoadAndKeep --readFilesIn ~/zjw/20190416upload/C1_CAGE-dir/trim_TATA_filter_yes/${i} --outFileNamePrefix ~/zjw/20190416upload/C1_CAGE-dir/mapping_output_filter_yes/${i}_ --outSAMtype SAM --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6
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


mkdir ~/zjw/20190416upload/C1_CAGE-dir/extract_uniquely_map_filter_yes/
for i in `ls ~/zjw/20190416upload/C1_CAGE-dir/mapping_output_filter_yes/ | grep "Aligned.out.sam"`
do
samtools view ~/zjw/20190416upload/C1_CAGE-dir/mapping_output_filter_yes/${i} | grep 'NH:i:1'$'\t''' > ~/zjw/20190416upload/C1_CAGE-dir/extract_uniquely_map_filter_yes/${i}_extract_uniquely_map.sam
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


mkdir ~/zjw/20190416upload/C1_CAGE-dir/split_plus_minus_filter_yes/
for i in `ls  ~/zjw/20190416upload/C1_CAGE-dir/extract_uniquely_map_filter_yes/ | grep "sam"`
do
cat ~/zjw/20190416upload/C1_CAGE-dir/extract_uniquely_map_filter_yes/${i} | awk '{FS=" "}{if ($2==0 || $2==256){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' > ~/zjw/20190416upload/C1_CAGE-dir/split_plus_minus_filter_yes/${i}_plus

cat ~/zjw/20190416upload/C1_CAGE-dir/extract_uniquely_map_filter_yes/${i} | awk '{FS=" "}{if ($2==16 || $2==272){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' > ~/zjw/20190416upload/C1_CAGE-dir/split_plus_minus_filter_yes/${i}_minus
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


mkdir ~/zjw/20190416upload/C1_CAGE-dir/extract_mismatch_filter_yes/

for i in `ls ~/zjw/20190416upload/C1_CAGE-dir/split_plus_minus_filter_yes/ | grep "extract_uniquely_map.sam_plus"`
do
python extractmismatch_plus.py -i ~/zjw/20190416upload/C1_CAGE-dir/split_plus_minus_filter_yes/${i} -o ~/zjw/20190416upload/C1_CAGE-dir/extract_mismatch_filter_yes/${i}_extractmismatch

python extractmismatch_minus.py -i ~/zjw/20190416upload/C1_CAGE-dir/split_plus_minus_filter_yes/${i%_*}_minus -o ~/zjw/20190416upload/C1_CAGE-dir/extract_mismatch_filter_yes/${i%_*}_minus_extractmismatch

cat ~/zjw/20190416upload/C1_CAGE-dir/extract_mismatch_filter_yes/${i}_extractmismatch    ~/zjw/20190416upload/C1_CAGE-dir/extract_mismatch_filter_yes/${i%_*}_minus_extractmismatch > ~/zjw/20190416upload/C1_CAGE-dir/extract_mismatch_filter_yes/${i%_*}_extractmismatch

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


mkdir ~/zjw/20190416upload/C1_CAGE-dir/final_out_filter_yes/
for i in `ls ~/zjw/20190416upload/C1_CAGE-dir/extract_mismatch_filter_yes/ | grep "sam_extractmismatch"`
do
samtools view -b -T ~/index/mm10_ERCC92/mm10_ERCC92.fa ~/zjw/20190416upload/C1_CAGE-dir/extract_mismatch_filter_yes/${i} | samtools view -b > ~/zjw/20190416upload/C1_CAGE-dir/final_out_filter_yes/${i}_add_header.bam

samtools sort ~/zjw/20190416upload/C1_CAGE-dir/final_out_filter_yes/${i}_add_header.bam -o ~/zjw/20190416upload/C1_CAGE-dir/final_out_filter_yes/${i}_add_header_sorted.bam

samtools index ~/zjw/20190416upload/C1_CAGE-dir/final_out_filter_yes/${i}_add_header_sorted.bam

bedtools bamtobed -i ~/zjw/20190416upload/C1_CAGE-dir/final_out_filter_yes/${i}_add_header_sorted.bam > ~/zjw/20190416upload/C1_CAGE-dir/final_out_filter_yes/${i}_add_header_sorted.bed
done
####################################################################################################################
