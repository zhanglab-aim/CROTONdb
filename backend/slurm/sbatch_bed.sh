#!/usr/bin/env bash
#SBATCH --time=2-1
#SBATCH --partition=gen

cd /mnt/ceph/users/vli/croton
export PYTHONPATH="."

python -u ./src/variant/make_beds.py --lstinx 7