#!/bin/bash

peak=$1
fa=$2
gtf_bed=$3
Rscript script/find_internal_TES_1.R ${peak}

python script/find_internal_TES.py ${fa} ${peak}_for_gs.bed > ${peak}_for_gs.bed.output

Rscript script/find_internal_TES_2.R ${peak} ${peak}_for_gs.bed.output ${gtf_bed}
