#!/bin/bash

mkdir ~/zjw/20190105
mkdir ~/zjw/20190105/3tail_read_with_tag
mkdir ~/zjw/20190105/3tail_read_with_tag_other_strand
mkdir ~/zjw/20190105/3tail_read_with_tag_other_strand_withA10_remain_A5
mkdir ~/zjw/20190105/mapping_output
mkdir ~/zjw/20190105/extract_uniquely_map
mkdir ~/zjw/20190105/split_plus_minus
mkdir ~/zjw/20190105/extract_mismatch
mkdir ~/zjw/20190105/add_header
mkdir ~/zjw/20190105/script_and_log

############### make sure that sample_list_tag.txt, extractmismatch_plus.py, extractmismatch_minus.py, and scCAT-seq_3'_data_processing.sh, cmpfastq_pe.pl, find_A10_and_remain_A5.py are in ~/zjw/20190105/script_and_log
############### make sure that fq files are list in ~/zjw/fastq_5cap_2018ab
############### STAR index must be built before running this script
############### reference genome fa file must be prepared in "add header" step

cd ~/zjw/20190105/script_and_log



for i in `ls ~/zjw/fastq_5cap_2018ab`
do
a=$(grep "${i}" ~/zjw/20190105/script_and_log/sample_list_tag.txt|awk '{print $2}')
if [[ ${a} =~ "NNNN" ]]; then
        b=${a#*4}
        echo ${b}
        echo "GTGGTATCAACGCAGAGT....${b%%TTTTTTTTTTTTTTTTTTTT*}"
        cat ~/zjw/fastq_5cap_2018ab/${i} | paste - - - - | grep "${b%%TTTTT*}TTTTT" | awk -v FS="\t" -v OFS="\n" '{print $1, $2, $3, $4}' > ~/zjw/20190105/3tail_read_with_tag/${i}_with_tag
fi
done



mkdir ~/zjw/20190914uploadCAT/tes_dir/other_strand_filter_yes/

for i in `ls ~/zjw/20190105/3tail_read_with_tag/ | grep "L1_1.fq_with_tag$"`
do
perl cmpfastq_pe.pl ~/zjw/20190105/3tail_read_with_tag/${i} ~/zjw/fastq_5cap_2018ab/${i%_1*}_2.fq
done

for i in `ls ~/zjw/20190105/3tail_read_with_tag/ | grep "L1_2.fq_with_tag$"`
do
perl cmpfastq_pe.pl ~/zjw/20190105/3tail_read_with_tag/${i} ~/zjw/fastq_5cap_2018ab/${i%_2*}_1.fq
done

rm ~/zjw/20190105/3tail_read_with_tag/*out
rm ~/zjw/fastq_5cap_2018ab/*unique.out
mv ~/zjw/fastq_5cap_2018ab/*out ~/zjw/20190105/3tail_read_with_tag_other_strand/






for i in `ls ~/zjw/20190105/3tail_read_with_tag_other_strand`
do
python ~/zjw/20190105/script_and_log/find_A10_and_remain_A5.py  -i ~/zjw/20190105/3tail_read_with_tag_other_strand/${i} -o ~/zjw/20190105/3tail_read_with_tag_other_strand_withA10_remain_A5/${i}_withA10_remain_A5
done



for i in `ls ~/zjw/20190105/3tail_read_with_tag_other_strand_withA10_trim/`
do
STAR --runThreadN 24 --genomeDir ~/index/mm10_ERCC92trimpolyA_STAR/ --genomeLoad LoadAndKeep --readFilesIn ~/zjw/20190105/3tail_read_with_tag_other_strand_withA10_trim/${i} --outFileNamePrefix ~/zjw/20190105/mapping_output/${i}_ --outSAMtype SAM --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6
done



for i in `ls ~/zjw/20190105/mapping_output/ |grep "sam"`
do
samtools view ~/zjw/20190105/mapping_output/${i} | grep 'NH:i:1'$'\t''' > ~/zjw/20190105/extract_uniquely_map/${i}_extract_uniquely_map.sam
done

for i in `ls ~/zjw/20190105/extract_uniquely_map | grep "sam"`
do
cat ~/zjw/20190105/extract_uniquely_map/${i} | awk '{FS=" "}{if ($2==0 || $2==256){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' > ~/zjw/20190105/split_plus_minus/${i}_plus
cat ~/zjw/20190105/extract_uniquely_map/${i} | awk '{FS=" "}{if ($2==16 || $2==272){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' > ~/zjw/20190105/split_plus_minus/${i}_minus
done


for i in `ls ~/zjw/20190105/split_plus_minus | grep "extract_uniquely_map.sam_plus"`
do
python ~/zjw/20190105/script_and_log/extractmismatch_plus.py -i ~/zjw/20190105/split_plus_minus/${i} -o ~/zjw/20190105/extract_mismatch/${i}_extractmismatch
python ~/zjw/20190105/script_and_log/extractmismatch_minus.py -i ~/zjw/20190105/split_plus_minus/${i%_*}_minus -o ~/zjw/20190105/extract_mismatch/${i%_*}_minus_extractmismatch
cat ~/zjw/20190105/extract_mismatch/${i}_extractmismatch ~/zjw/20190105/extract_mismatch/${i%_*}_minus_extractmismatch > ~/zjw/20190105/extract_mismatch/${i%_*}_extractmismatch
done

for i in `ls ~/zjw/20190105/extract_mismatch | grep "sam_extractmismatch"`
do
samtools view -b -T ~/index/mm10_ERCC92/mm10_ERCC92trimpolyA.fa ~/zjw/20190105/extract_mismatch/${i} | samtools view -b >  ~/zjw/20190105/add_header/${i}.bam
samtools sort ~/zjw/20190105/add_header/${i}.bam -o ~/zjw/20190105/add_header/${i}_sorted.bam
samtools index ~/zjw/20190105/add_header/${i}_sorted.bam
bedtools bamtobed -i ~/zjw/20190105/add_header/${i}_sorted.bam > ~/zjw/20190105/add_header/${i}.bed
done




for i in `ls  ~/zjw/20190105/add_header |grep "TKD"|grep "bed"|grep "L1_1"`
do
        a=$(wc -l ~/zjw/20190105/add_header/${i}|awk '{print $1}')
        b=$(wc -l ~/zjw/20190105/add_header/${i%%L1_1*}L1_2.fq-common.out_withA10_remain_A5_Aligned.out.sam_extract_uniquely_map.sam_extractmismatch.bed|awk '{print $1}')
        echo ${a}
        echo ${b}
        if [ ${a} -gt ${b} ]; then
                rm ~/zjw/20190105/add_header/${i%%L1_1*}L1_2*
        else
                rm ~/zjw/20190105/add_header/${i%%L1_1*}L1_1*
        fi
done
