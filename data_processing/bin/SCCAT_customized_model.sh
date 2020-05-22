#!/usr/bin/bash
set -e
FILE_PATH=$(cd "$(dirname "$0")"; pwd)
echo $FILE_PATH

ARGS=`getopt -a -o g:D:d:5:3:f:p:F:P:t:o:i:l: -l genome:,TrainTss:,TrainTes:,tss:,tes:,fantom5:,polyadb:,Fantom5:,Polyadb:,threshold:,output:,index:,library:,,help -- "$@"`
function usage() {
    echo "Usage:"
    echo ""
    echo "Description:"
    echo "  -g, --genome        Path to genome file."
    echo "  -D  --TrainTss      Directory of train TSS bed files."
    echo "  -d  --TrainTes      Directory of train TES bed files."
    echo "  -5, --tss           Directory of TSS bed files."
    echo "  -3, --tes           Directory of TES bed files."
    echo "  -f, --fantom5       Reference files from test FANTON5 databases."
    echo "  -F, --Fantom5       Reference files from train FANTON5 databases."
    echo "  -p, --polyadb       Reference files from test PolyA_DB3 databases."
    echo "  -P, --Polyadb       Reference files from train PolyA_DB3 databases."
    echo "  -t, --threshold     Threshold for filtering."
    echo "  -o, --output        Output filenames. Default = output"
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
      -D|--TrainTss)
              TrainTss="$2"
              shift;;
      -d|--TrainTes)
              TrainTes="$2"
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
      -F|--Fantom5)
              Fantom5="$2"
              shift;;
      -P|--Polyadb)
              Polyadb="$2"
              shift;;
      -t|--threshold)
              threshold="$2"
              shift;;
      -o|--output)
              output="$2"
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


echo genome:$genome:tss:$tss:tes$tes:fantom5:$fantom5:poly:$poly:threshold:$threshold:output:$output
echo TesModel:$TesModel:TssModel:$TssModel

if [ ! "${file##*.}"x = "bed"x ] && [ ! "${tes##*.}"x = "bed"x ]; then
        if [ $library = 10x ]; then
                sh script/fq2bed10xTSS.sh $genome $TrainTss $TrainTes $index
                sh script/fq2bed10xTES.sh $genome $TrainTss $TrainTes $index
        fi
        if [ $library = smart_seq ]; then
                sh script/fq2bedSCCATTSS.sh $genome $TrainTss $TrainTes $index
                sh script/fq2bedSCCATTES.sh $genome $TrainTss $TrainTes $index
        fi
        TrainTss=outdir/five_prime/collapse/sampleTSS_TKD.bed
        TrainTes=outdir/three_prime/collapse/sampleTES_TKD.bed
fi

sh $FILE_PATH/../script/callpeak_correction_TSS.sh $genome $TrainTss $Fantom5 $threshold
sh $FILE_PATH/../script/callpeak_correction_TES.sh $genome $TrainTes $Polyadb $threshold
wait # waiting for all tasks have benn finished

train_file=customized_model_training
python $FILE_PATH/../script/run_ml.py -e outdir/five_prime/peakfile/*.csv_final.csv -n $threshold -p $train_file
python $FILE_PATH/../script/run_ml.py -e outdir/three_prime/peakfile/*.csv_final.csv -n $threshold -p $train_file

TssModel=${train_file}*/tss/model/*.m
TesModel=${train_file}*/tes/model/*.m
rm -r outdir 

if [ -n "$tss" ]; then
        if [ -n "$tes" ]; then
                if [ ! "${file##*.}"x = "bed"x ] && [ ! "${tes##*.}"x = "bed"x ]; then
                        if [ $library = 10x ]; then
                                sh script/fq2bed10xTSS.sh $genome $tss $tes $index
                                sh script/fq2bed10xTES.sh $genome $tss $tes $index
                        fi
                        if [ $library = smart_seq ]; then
                                sh script/fq2bedSCCATTSS.sh $genome $tss $tes $index
                                sh script/fq2bedSCCATTES.sh $genome $tss $tes $index
                        fi
                        tss=outdir/five_prime/collapse/sampleTSS.bed
                        tes=outdir/three_prime/collapse/sampleTES.bed
                fi
                sh $FILE_PATH/../script/callpeak_correction_TSS.sh $genome $tss $fantom5 $threshold
                sh $FILE_PATH/../script/callpeak_correction_TES.sh $genome $tes $polyadb $threshold
		python $FILE_PATH/../script/run_ml.py -e outdir/five_prime/peakfile/*.csv_final.csv -m $TssModel -n $threshold
		python $FILE_PATH/../script/run_ml.py -e outdir/three_prime/peakfile/*.csv_final.csv -m $TesModel -n $threshold
                Rscript $FILE_PATH/../script/combine.R output_threshold*/tss/model/*_gene_rf_prediction.csv output_threshold*/tes/model/*_gene_rf_prediction.csv $output
                Rscript $FILE_PATH/../script/plot.R output_threshold*/ $output
                Rscript $FILE_PATH/../script/distance.R output_threshold*/tss/model/*_gene_rf_prediction.csv output_threshold*/tes/model/*_gene_rf_prediction.csv $output
                Rscript script/novel_temp.R output_threshold*/tss/model/*_gene_rf_prediction.csv output_threshold*/tss/model/*_intergenic_rf_prediction.csv output_threshold*/tes/model/*_gene_rf_prediction.csv output_threshold*/tes/model/*_intergenic_rf_prediction.csv $output
                Rscript $FILE_PATH/../script/pie.R output_threshold*/tss/model/*_gene_rf_prediction.csv output_threshold*/tss/model/*_intergenic_rf_prediction.csv output_threshold*/tes/model/*_gene_rf_prediction.csv output_threshold*/tes/model/*_intergenic_rf_prediction.csv $output

        fi
fi

# This part is used to predict the training data set itself when there is no other test set input
if [ -z "$tss" ]; then
        if [ -z "$tes" ]; then
                Rscript $FILE_PATH/../script/combine.R ${train_file}_threshold*/tss/model/*_gene_rf_prediction.csv ${train_file}_threshold*/tes/model/*_gene_rf_prediction.csv $output
                Rscript $FILE_PATH/../script/plot.R ${train_file}_threshold*/ $output
                Rscript $FILE_PATH/../script/distance.R ${train_file}_threshold*/tss/model/*_gene_rf_prediction.csv ${train_file}_threshold*/tes/model/*_gene_rf_prediction.csv $output
                Rscript script/novel_temp.R ${train_file}_threshold*/tss/model/*_gene_rf_prediction.csv ${train_file}_threshold*/tss/model/*_intergenic_rf_prediction.csv ${train_file}_threshold*/tes/model/*_gene_rf_prediction.csv ${train_file}_threshold*/tes/model/*_intergenic_rf_prediction.csv $output
                Rscript $FILE_PATH/../script/pie.R ${train_file}_threshold*/tss/model/*_gene_rf_prediction.csv ${train_file}_threshold*/tss/model/*_intergenic_rf_prediction.csv ${train_file}_threshold*/tes/model/*_gene_rf_prediction.csv ${train_file}_threshold*/tes/model/*_intergenic_rf_prediction.csv $output

        fi
fi

rm outdir output_threshold* -r