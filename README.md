
# scCAT-seq project  



Here is the code for our scCAT-seq project as described in the paper [scCAT-seq: single-cell identification and quantification of mRNA isoforms by cost-effective short-read sequencing of cap and tail](https://www.biorxiv.org/content/10.1101/2019.12.11.873505v1). scCAT-seq is designed to simultaneously capture transcript start sites(TSS) and transcript end sites(TES) on single-cell level.

---

# The pipeline contains two parts:
* The first part contains the scripts for data pre-processing, TSS/TES peak calling, and correction by machine learning algorithms. The scripts and demos are included in the ["data_processing" directory](https://github.com/huyoujinlab/scCAT-seq/tree/master/data_processing).

* The second part contains the scripts to further analyze the data and generate the relevant figures, which are deposited in the [result directory](https://github.com/huyoujinlab/scCAT-seq/tree/master/result)
