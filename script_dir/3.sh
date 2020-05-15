#!/bin/bash

peak=$1
fa=$2
gtf_bed=$3
Rscript 3_1.R ${peak}

python 3.py ${fa} ${peak}_for_gs.bed > ${peak}_for_gs.bed.output

Rscript 3_2.R ${peak} ${peak}_for_gs.bed.output ${gtf_bed}
