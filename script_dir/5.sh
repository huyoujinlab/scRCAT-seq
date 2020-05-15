#!/bin/bash

set -e

peak=$1
fa=$2
gtf_bed=$3

Rscript 5_1.R ${peak}

python 5.py ${fa} ${peak}_for_gs.bed > ${peak}_for_gs.bed.output

Rscript 5_2.R ${peak} ${peak}_for_gs.bed.output ${gtf_bed}
