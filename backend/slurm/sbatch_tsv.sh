#!/usr/bin/env bash
#SBATCH --time=2-1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100-32gb:1

cd /mnt/ceph/users/vli/croton
export PYTHONPATH="."

python -u ./src/variant/make_tsvs.py --lstinx 7