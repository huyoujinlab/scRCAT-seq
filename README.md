
# scCAT-seq project  



Here is the code for our scCAT-seq project as described in the paper [scCAT-seq: single-cell identification and quantification of mRNA isoforms by cost-effective short-read sequencing of cap and tail](https://www.biorxiv.org/content/10.1101/2019.12.11.873505v1). scCAT-seq is designed to simultaneously capture transcript start sites(TSS) and transcript end sites(TES) on single-cell level.


### This project consists of two parts.

* The first part contains the scripts for data pre-processing, TSS/TES peak calling, and correction by machine learning algorithms. The scripts and demos are included in the ["data_processing" directory](https://github.com/huyoujinlab/scCAT-seq/tree/master/data_processing).

* The second part contains the scripts to further analyze the data and generate the relevant figures, which are deposited in the ["result" directory] (https://github.com/huyoujinlab/scCAT-seq/tree/master/result).

Readers can switch to the subfolder to get detailed information.

---
### The softwares and packages required:

On linux:


1) Python 2.7.15   
2) R 3.6.1  
3) Perl 5.26.2
3) STAR 2.7.3a   
4) Samtools 1.3.1   
5) Bedtools 2.27.1   
6) cDNA_cupcake 8.0   
7) Cutadapt 1.18   
8) Minimap2 2.17   
9) Cell ranger 3.1.0   
 
 
---

On R 3.6.1:
 
1) basicTrendline 2.0.3   
2) broom 0.5.2
3) BSgenome 1.54.0
4) BSgenome.Ercc92.ZJW.1 1.0
5) BSgenome.Mmusculus.UCSC.mm10 1.4.0
6) BSgenome.Hsapiens.UCSC.hg38 1.4.1
7) BuenColors 0.5.5   
8) CAGEr 1.24.0   
9) data.table 1.12.4   
10) DESeq2 1.22.1   
11) dplyr 0.8.3
12) genefilter 1.64.0
13) ggalt 0.4.0   
14) ggplot2 3.2.1   
15) ggpubr 0.2.3   
16) ggsignif 0.6.0   
17) gplots 3.0.1.1   
18) gridExtra 2.3   
19) heatmap.plus 1.3   
20) Monocle2 2.12.0   
21) nlstools 1.0-2   
22) pheatmap 1.0.12   
23) purrr 0.3.2   
24) RColorBrewer 1.1-2   
25) reshape 0.8.8   
26) rlist 0.4.6.1   
27) Rmisc 1.5   
28) rsample 0.0.5   
29) rstatix 0.2.0   
30) scales 1.0.0   
31) splines 3.5.1 
32) statmod 1.4.32
33) stringr 1.4.3   
34) Seurat 3.0.2   
 
BS.genome packages were needed for `CAGEr`. As BS.genome packages for ERCC are not available publicly, we built BS.genome package for ERCC and uploaded at (https://drive.google.com/open?id=1cwJSUWcZ8PkYAs7vmUGDYpI_fW3jE3nK). 

---
On python 3.7:
 
1) pandas 0.24.1  
2) regex 2019.6.5
3) numpy 1.16.5  
4) matplotlib 2.2.3  
5) seaborn 0.8.1
6) sklearn 0.21.3  
7) scikit-plot 0.3.7  


---


### To install the softwares and packages:


```
### Download anaconda
wget https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh

### Install anaconda
sh Anaconda2-2019.10-Linux-x86_64.sh


1) Creating environment and installing software
### Bulid environment
conda create -n scCAT_seq python=3.7
conda activate scCAT_seq

### Install softwares and packages
conda install -c anaconda perl=5.26.2
conda install -c bioconda STAR=2.7.3a
conda install -c bioconda samtools=1.3.1
conda install -c bioconda bedtools=2.27.1
conda install -c bioconda cutadapt=1.18

conda install -c conda-forge pandas=0.24.1
conda install -c conda-forge regex
conda install -c conda-forge plotnine
conda install -c conda-forge python-levenshtein
conda install -c conda-forge pandas=0.24.1
conda install -c conda-forge scikit-learn=0.22.1

conda install -c r r-base=3.6.1
conda install -c conda-forge r-foreign
conda install -c conda-forge r-ggal
Rscript install_R_packages.R


2) Or readers can creat environment from scCAT_seq.yml file
conda env create -f scCAT_seq.yml
Rscript install_R_packages.R
```

