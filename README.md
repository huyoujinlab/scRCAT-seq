
# scRCAT-seq project  



* Here is the code for our scRCAT-seq project as described in the paper [scCAT-seq: single-cell identification and quantification of mRNA isoforms by cost-effective short-read sequencing of cap and tail](https://www.biorxiv.org/content/10.1101/2019.12.11.873505v1). scRCAT-seq is designed to simultaneously capture RNA transcript start sites(TSSs) and transcript end sites(TESs) on single-cell level.

---

# The pipeline contains two parts:
* The first part contains the scripts for data pre-processing, TSS/TES peak calling, and correction by machine learning algorithms. The scripts and demos are included in the [`data_processing`](https://github.com/huyoujinlab/scRCAT-seq/tree/master/data_processing) directory.

* The second part contains the scripts to further analyze the data and generate the relevant figures, which are deposited in the [`result` ](https://github.com/huyoujinlab/scRCAT-seq/tree/master/result) directory.
