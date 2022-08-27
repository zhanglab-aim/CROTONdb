#!/usr/bin/env bash
#SBATCH --time=2-1
#SBATCH --partition=ccb

cd /mnt/ceph/users/vli/croton
export PYTHONPATH="."

python -u ./src/variant/make_tabix.py --lstinx 7