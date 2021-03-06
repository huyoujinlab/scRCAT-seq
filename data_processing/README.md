# Description

The pipeline is used to process the scRCAT-seq data to identify TSSs/TESs of transcripts based on random forests (RF) model. With the model trained on the dataset derived from human embryonic stem cells, for which TSSs and TESs are well annotated in the `FANTOM5` database and `PolyA_DB` respectively, the pipeline can be used to identify authentic transcript boundaries (TSSs and TESs) for cells and species with scRCAT-seq data. Alternatively, users can use their own training datasets to train the model, which can then be used to identify TSSs/TESs from their own dataset.
The pipeline can calculate the accuracy to identify authentic TSSs/TESs, and identify novel TSSs and TESs, novel genes for single cells. The number of alternative TSSs/TESs, and major TSS and major TES of each gene can also be identified.

# System requirements

The pipeline has been tested on centos and Ubuntu operating systems, and the software required are listed below.

## Installation


Users can install softwares and packages by conda

### Download anaconda

```
wget https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh
```

### Install anaconda

```
sh Anaconda2-2019.10-Linux-x86_64.sh
```


The softwares can be installed in two ways:

1) Use the following commands:
   
```
### Build environment
conda create -n scCAT_seq python=3.7
conda activate scCAT_seq

### Install the softwares and packages
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
conda install -c conda-forge r-base=3.6.1
conda install -c conda-forge r-foreign
conda install -c conda-forge r-ggal
Rscript install/install_R_packages.R
```

2) Or readers can creat environment from `scCAT_seq.yml` file (recommended)

```
conda env create -f install/scRCAT_seq.yml
conda activate scRCAT_seq
Rscript install/install_R_packages.R
```
   
This step will spend 18 minutes to install all software and packages.



# Demo

Here, we demo the pipeline to identify TSSs and TESs for human embryonic stem cells. 

## Demo data download

```
sh script/download.sh
```

## Usage:

```
sh ./bin/SCCAT.sh -g ~/index/hg38/hg38.fa -5 input/hESCnofiltersccatUMI5_TKD.bed -3 input/hESCnofiltersccatUMI3_TKD.bed -f reference/tc_hESC.bed -p reference/human.PAS100_hg38.bed -t 3 -o output
```

Parameters in the pipeline are described as follows: 

Input files for this pipeline are the bed files of TESs and TSSs, which was generated by mapping fastq files with `STAR` and transforming with `samtools`.

`-g`  The filename of genome, which is used to find motif around TESs and TESs

`-5` Path to the bed file containing TSS peaks or R1 paired-end fastq files.

`-3` Path to the bed file containing TES peaks or R2 paired-end fastq files.

`-f` Reference to annotate authentic TSSs of the data to be predicted (specified by -5).

`-p` Reference to annotate authentic TESs of the data to be predicted (specified by -3).

`-t` Threshold for filtering out the peaks of low abundance.

`-o` Folder for storing output files.

Use `SCCAT.sh –h` to get more detailed usage. 

This step will take 16min to finish.

## Output:

In this demo pipeline, output files were deposited in the `output/` directory with default. Including `tc_hESCnofiltersccatUMI3_5cap_prediction.tsv`, `tc_hESCnofiltersccatUMI3_3tail_prediction.tsv`, `combine_result.tsv`, `novelgene.tsv`, `accuracy.pdf`, `distance.pdf`, `pieTSS.pdf`, `pieTES.pdf`, `novel.pdf`.

The first file is `tc_hESCnofiltersccatUMI5_5cap_prediction.tsv`. The predicted authentic TSSs from the candidate peaks. Each row represents a peak. The colnames represent:

1)	gene: gene symbol of the peak
2)	chr: chromosome information of the peak
3)	peak_position: the position of the peak
4)	strand: strand of the peak
5)	model_prediction: Predicted value, the value 1 indicates a true TSS peak while value 0 indicates a false TSS peak.


The second file is `tc_hESCnofiltersccatUMI3_3tail_prediction.tsv`. The predicted authentic TESs from the candidate peaks. Each row represents a peak. The colnames represent:

1)	gene: gene symbol of the peak
2)	chr: chromosome information of the peak
3)	peak_pos: the position of the peak
4)	strand: strand information of the peak
5)	model.prediction: Predicted value, the value 1 indicates a true TES peak while value 0 indicates a false TES peak.

The third file is `combine_result.tsv`. TSS/TES choices and dominant TSSs/TESs of each gene. Each row represents a gene. The colnames represent:

1)	gene: gene symbol name
2)	chr: chromosome information of the gene
3)	strand: strand information of the gene
4)	Num_TSS: Number of authentic TSSs identified for the gene.
5)	Num_TES: Number of authentic TESs identified for the gene.
6)	TSS_pos: Position of the dominant TSS.
7)	TES_pos: Position of the dominant TES.


The forth file is `novelgene.tsv`. Each row represents a novel gene. Each column represent:
 
1)	chr: chromosome information of the novel gene
2)	start: start position of the novel gene
3)	end: end position of the novel gene
4)	strand: strand information of novel gene

Plots:

The fifth file is `accuracy.pdf`. It shows the accuracy in identifying authentic TSSs and TESs with rf machine learning models.
 
The sixth file is `distance.pdf`. It shows the distance of the identified TSSs/TESs to those annotated in hg38. 
 
The seventh file is `pieTSS.pdf`. Pie chart shows the genomic distribution of the identified TSSs in hESC. 
 
The eighth file is `pieTES.pdf`. Pie chart shows the genomic distribution of the identified TESs in hESC. 

The ninth file is `novel.pdf`. It shows the number of novel isoforms of annotated genes and novel, unannotated transcripts in hESC.

 

# Training on the model with own data (optional)
In some cases, users may use their own dataset to train the model, they can use the `SCCAT_customized_model.sh`, when dataset for training should be deposit in the `customized_model_training_threshold3` file. The model will be trained with 70% of the data and tested on the rest 30%. The accuracy for the testing data would be saved in the `customized_model_training_threshold3/` folder. 
In the meantime, other dataset need to be processed to identify TSS and TES are input by specifying the directory after the paremeter `-5` and  `-3` respectively. The output files would be found in the output files.

For example, here we train the model with hESC dataset, and use the model to predict the TSS and TES for same dataset. Users can choose different datasets. 

```
sh ./bin/SCCAT_customized_model.sh -g ~/index/hg38/hg38.fa  -D input/hESCnofiltersccatUMI5_TKD.bed -d input/hESCnofiltersccatUMI3_TKD.bed -5 input/hESCnofiltersccatUMI5_TKD.bed -3 input/hESCnofiltersccatUMI3_TKD.bed -F reference/tc_hESC.bed -P reference/human.PAS100_hg38.bed -f reference/tc_hESC.bed -p reference/human.PAS100_hg38.bed -t 3 -o output_customized
```

`-g` The filename of the reference genome

`-D` Path to the TSS training data.

`-d` Path to the TES training data.

`-5` Directory of the TSS data to be predicted.

`-3` Directory of the TES data to be predicted.

`-F` Reference to annotate authentic TSSs of training data.

`-P` Reference to annotate authentic TESs of training data.

`-f` Reference to annotate authentic TSSs of the data to be predicted (specified by -5).

`-p` Reference to annotate authentic TESs of the data to be predicted (specified by -3).

`-t` Threshold for filtering out peaks with low abundance.

`-o` Directory to save output files.

Output:

1. Accuracy for the customized models would be found in the `customized_model_training_threshold3/` folder.
2. The following results would be saved in the directory `output_customized/`: the accuracy to identify authentic TSSs/TESs, and novel TSSs and TESs, novel genes for single cells. The number of alternative TSSs/TESs, and major TSS and TES of each gene also be identified.
