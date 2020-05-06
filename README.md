
# scCAT-seq project  



Here is the code for our scCAT-seq project as described in the paper [scCAT-seq: single-cell identification and quantification of mRNA isoforms by cost-effective short-read sequencing of cap and tail](https://www.biorxiv.org/content/10.1101/2019.12.11.873505v1). scCAT-seq is designed to simultaneously capture transcript start sites(TSS) and transcript end sites(TES) on single-cell level.


### This project consists of two parts.

* The [first part](https://github.com/huyoujinlab/scCAT-seq/tree/master/data_processing) is code used to data pre-processing and peak correction.

* The [second part](https://github.com/huyoujinlab/scCAT-seq/tree/master/Analysing) is code used to analyse the data from scCAT-seq.

You can switch to the subfolder to get detailed information.

---
### Some softwares and packages used in this publication are listed below.

On linux:


1) Python 2.7.15   
2) R 3.5.0   
3) STAR 2.6.1a   
4) Samtools 1.3.1   
5) Bedtools 2.27.1   
6) cDNA_cupcake 8.0   
7) Cutadapt 1.18   
8) Minimap2 2.17   
9) Cell ranger 2.1.1   
 
 
---

On R 3.5.0:
 
1) basicTrendline 2.0.3   
2) broom 0.5.2
3) BSgenome 1.54.0
4) BSgenome.Ercc92.ZJW.1 1.0
5) BSgenome.Macfas.ZJW.1 1.0
6) BSgenome.Mmusculus.UCSC.mm10 1.4.0
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
 
BS.genome packages were needed for `CAGEr`. However, BS.genome packages for ERCC and crab-eating monkey are no available. We have built the packages. User can download and install for reproducting out results.

* BS.genome package for ERCC is provided [here](https://drive.google.com/open?id=1cwJSUWcZ8PkYAs7vmUGDYpI_fW3jE3nK).

* BS.genome package for crab-eating monkey is provided [here](https://drive.google.com/open?id=1ZvdSOCV5AbwcqYPZ_SYX42egSRKUGBrC).
---
On python 2.7.15:
 
1) pandas 0.24.1  
2) regex 2019.6.5
3) numpy 1.16.5  
4) matplotlib 2.2.3  
5) seaborn 0.8.1
6) sklearn 0.21.3  
7) scikit-plot 0.3.7  


---


Reader can install softwares and packages by conda


```
### Download anaconda
wget https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh

### Install anaconda
sh Anaconda2-2019.10-Linux-x86_64.sh

### Bulid environment
conda create -n scCAT_seq python=2.7.15
conda activate scCAT_seq

### Install softwares and packages
conda install -c bioconda STAR=2.6.1a
conda install -c bioconda samtools=1.3.1
conda install -c bioconda bedtools=2.27.1
conda install -c bioconda cutadapt=1.18
conda install -c conda-forge r-base 

conda install -c conda-forge pandas=0.24.1
conda install -c conda-forge regex=2019.6.5
conda install -c conda-forge numpy=1.16.5
conda install -c conda-forge matplotlib=2.2.3
conda install -c conda-forge seaborn=0.8.1
conda install -c conda-forge scikit-plot=0.3.7
```
