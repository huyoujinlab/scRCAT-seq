# Demo for machine learning


#### Here is a demo to call peak and machine learning. After we generate bed file by scCAT-seq_5_data_processing.sh/scCAT-seq_3_data_processing.sh, we use them to call peak. We extract the reads on chr 3 for simplification



#### First, we run:

```
bash CAGE_dominant.sh 
```

* "D44_52_3tail_dominant_tes.bed" and "D44_52_5cap_dominant_tss.bed" file will be generated. These file is peak calling output.



#### Second, we run:

```
Rscript cal_slope_cattss.R D44_52_5cap_dominant_tss_upstream2k_and_genebody.bed D_merge_standard_smart_seq2.count
Rscript cal_slope_cattes.R D44_52_3tail_dominant_tes_genebody_downstream2k.bed D_merge_standard_smart_seq2.count
```

* "tc_D44_52_5cap_peak.csv" and "tc_D44_52_3tail_peak.csv" will be generated. In this sept, we added some features for machine learning.



#### Third, we run:

```
Rscript cal_slope_cattss2.R tc_D44_52_5cap_peak.csv D_merge_standard_smart_seq2_sorted_chr3.wig
Rscript cal_slope_cattes2.R tc_D44_52_3tail_peak.csv  D_merge_standard_smart_seq2_sorted_chr3.wig
```

* "tc_D44_52_5cap_peak.csvnew.csv" and "tc_D44_52_3tail_peak.csvnew.csv" will be generated. In this sept, we calculate the slope of Smart-seq2 coverage curve around the peaks.



#### Fourth, we run:

```
Rscript cal_slope_cattss3.R tc_D44_52_5cap_peak.csvnew.csv
Rscript cal_slope_cattes3.R tc_D44_52_3tail_peak.csvnew.csv
```

* "tc_D44_52_5cap_peak.csvnew.csv.csv" and "tc_D44_52_3tail_peak.csvnew.csv.csv" will be generated. In this sept, features were added and we can define the TRUE or FALSE of each peak.





#### Fifth, we run:

```
mv tc_D44_52_5cap_peak.csvnew.csv.csv ../data/
mv tc_D44_52_3tail_peak.csvnew.csv.csv ../data/
cd ../
python ./bin/demo.py
```

* The new folder named "output" will be built after completion of the pipeline, which contains all the output files, including the "\*.csv" file. 

* The last column in the .CSV files shows the predicted classification corresponding to each individual peaks in the input file. The value "1" indicates a true TSS/TES peak while value "0" indicates a false TSS/TES peak.


* from sept 1 to sept 5, it takes about 10 minutes.







