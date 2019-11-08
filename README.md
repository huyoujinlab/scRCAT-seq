Code for our scCAT-seq project



Some softwares used in this study on linux paltform was listed below:

Software on linux:
Python 2.7.15   
R 3.5.0   
STAR 2.6.1a   
Samtools 1.3.1   
Bedtools 2.27.1   
cDNA_cupcake 8.0   
Cutadapt 1.18   
Minimap2 2.17   
Cell ranger 2.1.1   

 
Some packages used in this study on R was listed below:

basicTrendline 2.0.3   
broom 0.5.2   
BuenColors 0.5.5   
CAGEr 1.24.0   
data.table 1.12.4   
DESeq2 1.22.1   
dplyr 0.8.3   
ggalt 0.4.0   
ggplot2 3.2.1   
ggpubr 0.2.3   
ggsignif 0.6.0   
gplots 3.0.1.1   
gridExtra 2.3   
heatmap.plus 1.3   
Monocle2 2.12.0   
nlstools 1.0-2   
pheatmap 1.0.12   
purrr 0.3.2   
RColorBrewer 1.1-2   
reshape 0.8.8   
rlist 0.4.6.1   
Rmisc 1.5   
rsample 0.0.5   
rstatix 0.2.0   
scales 1.0.0   
splines 3.5.1   
stringr 1.4.3   
Seurat 3.0.2   




Here is a demo to call peak and machine learning

### after we generate bed file by scCAT-seq_5_data_processing.sh/scCAT-seq_3_data_processing.sh, we use them to call peak.

### We extract the reads on chr 3 for simplification



### First, we run:

bash CAGE_dominant.sh 

### "D44_52_3tail_dominant_tes.bed" and "D44_52_5cap_dominant_tss.bed" file will be generated. These file is peak calling output.



### Second, we run:

Rscript cal_slope_cattss.R D44_52_5cap_dominant_tss_upstream2k_and_genebody.bed D_merge_standard_smart_seq2.count
Rscript cal_slope_cattes.R D44_52_3tail_dominant_tes_genebody_downstream2k.bed D_merge_standard_smart_seq2.count

### "tc_D44_52_5cap_peak.csv" and "tc_D44_52_3tail_peak.csv" will be generated. In this sept, we added some features for machine learning.



### Third, we run:

Rscript cal_slope_cattss2.R tc_D44_52_5cap_peak.csv D_merge_standard_smart_seq2_sorted_chr3.wig
Rscript cal_slope_cattes2.R tc_D44_52_3tail_peak.csv  D_merge_standard_smart_seq2_sorted_chr3.wig

### "tc_D44_52_5cap_peak.csvnew.csv" and "tc_D44_52_3tail_peak.csvnew.csv" will be generated. In this sept, we calculate the slope of Smart-seq2 coverage curve around the peaks.



### Fourth, we run:

Rscript cal_slope_cattss3.R tc_D44_52_5cap_peak.csvnew.csv
Rscript cal_slope_cattes3.R tc_D44_52_3tail_peak.csvnew.csv

### "tc_D44_52_5cap_peak.csvnew.csv.csv" and "tc_D44_52_3tail_peak.csvnew.csv.csv" will be generated. In this sept, features were added and we can define the TRUE or FALSE of each peak.




WARNING: Initial version of code, still a work in progress.
