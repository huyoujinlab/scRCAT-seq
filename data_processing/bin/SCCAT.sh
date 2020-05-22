#!/usr/bin/bash
set -e
FILE_PATH=$(cd "$(dirname "$0")"; pwd)
echo $FILE_PATH

ARGS=`getopt -a -o g:5:3:f:p:t:o:m:M:i:l: -l genome:,tss:,tes:,fantom5:,polyadb:,threshold:,output:,TssModel:,TesModel:,index:,library:,help -- "$@"`
function usage() {
    echo "Usage:"
    echo ""
    echo "Description:"
    echo "  -g, --genome        Path to genome file, which is used to find motif around TESs and TESs."
    echo "  -5, --tss           Path to the bed file containing TSS peaks or R1 paired-end fastq files."
    echo "  -3, --tes           Path to the bed file containing TES peaks or R2 paired-end fastq files."
    echo "  -f, --fantom5       Reference to annotate authentic TSSs of the data to be predicted (specified by -5)."
    echo "  -p, --polyadb       Reference to annotate authentic TESs of the data to be predicted (specified by -3)."
    echo "  -t, --threshold     Threshold for filtering out the peaks of low abundance."
    echo "  -o, --output        Folder for storing output files. default = output"
    echo "  -m, --TssModel      Path to tss model file."
    echo "  -M, --TesModel      Path to tes model file."
    echo "  -i, --index         Path to STAR index folder. This parameter is required when the input data is a fastq file."
    echo "  -l, --library       How to build a library, we have two method: 10x or smart_seq, defacult = smart_seq. This parameter is required when the input data is a fastq file."
    exit -1   
}
[ $? -ne 0 ] && usage
#set -- "${ARGS}"
eval set -- "${ARGS}"
while true
do
      case "$1" in
      -g|--genome)
              genome="$2"
              shift;;
      -5|--tss)
              tss="$2"
              shift;;
      -3|--tes)
              tes="$2"
              shift;;
      -f|--fantom5)
              fantom5="$2"
              shift;;
      -p|--polyadb)
              polyadb="$2"
              shift;;
      -t|--threshold)
              threshold="$2"
              shift;;
      -o|--output)
              output="$2"
              shift;;
      -m|--TssModel)
              TssModel="$2"
              shift;;
      -M|--TesModel)
              TesModel="$2"
              shift;; 
      -i|--index)
              index="$2"
              shift;; 
      -l|--library)
              library="$2"
              shift;; 
      -h|--help)
              usage;;
      --)
              shift
              break;;
      esac
shift
done 


if [ -z "$output" ]; then
    output=output
fi

if [ -z "$threshold" ]; then
    threshold=3
fi

if [ -z "$library" ]; then
    library=smart_seq
fi

if [ -z "$TssModel" ]; then
    if [ $library = 10x ]; then
        TssModel=$FILE_PATH/../reference/tc_organiodnofilter10xUMI5threshold0_5cap_final_gene_rf_model.m
    fi
    if [ $library = smart_seq ]; then
        TssModel=$FILE_PATH/../reference/hESC_training_rf_model_tss.m
    fi
fi

if [ -z "$TesModel" ]; then
    if [ $library = 10x ]; then
        TesModel=$FILE_PATH/../reference/tc_D80_3_3tail_final_gene_rf_model.m
    fi
    if [ $library = smart_seq ]; then
        TesModel=$FILE_PATH/../reference/hESC_training_rf_model_tes.m
    fi
fi

echo TesModel:$TesModel:TssModel:$TssModel

if [ ! "${tss##*.}"x = "bed"x ] && [ ! "${tes##*.}"x = "bed"x ]; then
        if [ $library = 10x ]; then
                sh script/fq2bed10xTSS.sh $genome $tss $tes $index
                sh script/fq2bed10xTES.sh $genome $tss $tes $index
        fi
        if [ $library = smart_seq ]; then
                sh script/fq2bedSCCATTSS.sh $genome $tss $tes $index
                sh script/fq2bedSCCATTES.sh $genome $tss $tes $index
        fi
        tss=outdir/five_prime/collapse/sampleTSS_TKD.bed
        tes=outdir/three_prime/collapse/sampleTES_TKD.bed
fi


sh $FILE_PATH/../script/callpeak_correction_TSS.sh $genome $tss $fantom5 $threshold
sh $FILE_PATH/../script/callpeak_correction_TES.sh $genome $tes $polyadb $threshold

wait # waiting for all tasks have benn finished

python $FILE_PATH/../script/run_ml.py -e outdir/five_prime/peakfile/*.csv_final.csv -m $TssModel -n $threshold
python $FILE_PATH/../script/run_ml.py -e outdir/three_prime/peakfile/*.csv_final.csv -m $TesModel -n $threshold


Rscript $FILE_PATH/../script/combine.R output_threshold*/tss/model/*_gene_rf_prediction.csv output_threshold*/tes/model/*_gene_rf_prediction.csv $output
Rscript $FILE_PATH/../script/plot.R output_threshold*/ $output
Rscript $FILE_PATH/../script/distance.R output_threshold*/tss/model/*_gene_rf_prediction.csv output_threshold*/tes/model/*_gene_rf_prediction.csv $output
Rscript script/novel_temp.R output_threshold*/tss/model/*_gene_rf_prediction.csv output_threshold*/tss/model/*_intergenic_rf_prediction.csv output_threshold*/tes/model/*_gene_rf_prediction.csv output_threshold*/tes/model/*_intergenic_rf_prediction.csv $output
Rscript $FILE_PATH/../script/pie.R output_threshold*/tss/model/*_gene_rf_prediction.csv output_threshold*/tss/model/*_intergenic_rf_prediction.csv output_threshold*/tes/model/*_gene_rf_prediction.csv output_threshold*/tes/model/*_intergenic_rf_prediction.csv $output
rm outdir *threshold* -r