#!/usr/bin/env bash
#SBATCH --time=2-1
#SBATCH --partition=defq

cd /home/lingj/zhanglab/jling/CROTONdb/backend
#export PYTHONPATH="."

for i in `seq 1 12`; do echo $i; python make_beds.py --chrom $i; done 

#python -u make_beds.py --lstinx 7