#!/bin/bash

mkdir ~/zjw/20190109
mkdir ~/zjw/20190109/5cap_read_with_tag
mkdir ~/zjw/20190109/trim_GTGGTATCAACGCAGAGTACAT
mkdir ~/zjw/20190109/mapping_output
mkdir ~/zjw/20190109/extract_uniquely_map
mkdir ~/zjw/20190109/split_plus_minus
mkdir ~/zjw/20190109/extract_mismatch
mkdir ~/zjw/20190109/final_out
mkdir ~/zjw/20190109/script_and_log

############### make sure that extractmismatch_plus.py, extractmismatch_minus.py, and scCAT-seq_5'_data_processing.sh are in ~/zjw/20190109/script_and_log
############### make sure that fq files are list in ~/zjw/fastq_5cap_2018ab
############### STAR index must be built before running this script
############### reference genome fa file must be prepared in "add header" step

cd ~/zjw/20190109/script_and_log

cd ~/zjw/20190109/script_and_log

echo "###############read with tag"
for i in `ls ~/zjw/fastq_5cap_2018ab/`
do
a=$(grep "${i}" ~/zjw/20190109/script_and_log/sample_list_tag.txt|awk '{print $6}')
echo ${a}
cat ~/zjw/fastq_5cap_2018ab/${i} | paste - - - - | grep $'\t'"${a}" | awk -v FS="\t" -v OFS="\n" '{print $1, $2, $3, $4}' > ~/zjw/20190109/5cap_read_with_tag/${i}_with_tag.fq
done

echo "###############trim GTGGTATCAACGCAGAGTACAT"
for i in `ls ~/zjw/20190109/5cap_read_with_tag`
do
cutadapt -g GTGGTATCAACGCAGAGTACAT -o ~/zjw/20190109/trim_GTGGTATCAACGCAGAGTACAT/${i}.trimed.remainGGG  ~/zjw/20190109/5cap_read_with_tag/${i}
done

echo "###############mapping"
for i in `ls ~/zjw/20190109/trim_GTGGTATCAACGCAGAGTACAT/`
do
STAR --runThreadN 24 --genomeDir ~/index/mm10_ERCC92trimpolyA_STAR/ --genomeLoad LoadAndKeep --readFilesIn ~/zjw/20190109/trim_GTGGTATCAACGCAGAGTACAT/${i} --outFileNamePrefix ~/zjw/20190109/mapping_output/${i}_ --outSAMtype SAM --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6
done

echo "###############extract uniquely map"
for i in `ls ~/zjw/20190109/mapping_output|grep "sam"`
do
samtools view ~/zjw/20190109/mapping_output/${i} |grep 'NH:i:1'$'\t''' > ~/zjw/20190109/extract_uniquely_map/${i}_extract_uniquely_map.sam
done

echo "###############split plus minus"
for i in `ls ~/zjw/20190109/extract_uniquely_map | grep "sam"`
do
cat ~/zjw/20190109/extract_uniquely_map/${i} | awk '{FS=" "}{if ($2==0 || $2==256){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' > ~/zjw/20190109/split_plus_minus/${i}_plus
cat ~/zjw/20190109/extract_uniquely_map/${i} | awk '{FS=" "}{if ($2==16 || $2==272){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' > ~/zjw/20190109/split_plus_minus/${i}_minus
done

echo "###############extract mismatch"
for i in `ls ~/zjw/20190109/split_plus_minus | grep "extract_uniquely_map.sam_plus"`
do
python ~/zjw/20190109/script_and_log/extractmismatch_plus.py -i ~/zjw/20190109/split_plus_minus/${i} -o ~/zjw/20190109/extract_mismatch/${i}_extractmismatch
python ~/zjw/20190109/script_and_log/extractmismatch_minus.py -i ~/zjw/20190109/split_plus_minus/${i%_*}_minus -o ~/zjw/20190109/extract_mismatch/${i%_*}_minus_extractmismatch
cat ~/zjw/20190109/extract_mismatch/${i}_extractmismatch ~/zjw/20190109/extract_mismatch/${i%_*}_minus_extractmismatch > ~/zjw/20190109/extract_mismatch/${i%_*}_extractmismatch
done

echo "###############final out"
for i in `ls ~/zjw/20190109/extract_mismatch | grep "sam_extractmismatch"`
do
samtools view -b -T ~/index/mm10_ERCC92/mm10_ERCC92trimpolyA.fa ~/zjw/20190109/extract_mismatch/${i} | samtools view -b >  ~/zjw/20190109/final_out/${i}_add_header.bam
samtools sort ~/zjw/20190109/final_out/${i}_add_header.bam -o ~/zjw/20190109/final_out/${i}_add_header_sorted.bam
samtools index ~/zjw/20190109/final_out/${i}_add_header_sorted.bam
bedtools bamtobed -i ~/zjw/20190109/final_out/${i}_add_header_sorted.bam > ~/zjw/20190109/final_out/${i}_add_header_sorted.bed
done

echo "###############remove one end"
for i in `ls  ~/zjw/20190109/final_out |grep "TKD"|grep "bed"|grep "L1_1"`
do
        a=$(wc -l ~/zjw/20190109/final_out/${i}|awk '{print $1}')
        b=$(wc -l ~/zjw/20190109/final_out/${i%%L1_1*}L1_2.fq_with_tag.fq.trimed.remainGGG_Aligned.out.sam_extract_uniquely_map.sam_extractmismatch_add_header_sorted.bed|awk '{print $1}')
        echo ${a}
        echo ${b}
        if [ ${a} -gt ${b} ]; then
                rm ~/zjw/20190109/final_out/${i%%L1_1*}L1_2*
        else
                rm ~/zjw/20190109/final_out/${i%%L1_1*}L1_1*
        fi
done
