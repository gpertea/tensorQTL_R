#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=64G,h_vmem=64G,h_fsize=100G
#$ -N tqtlMDD_Amyg
#$ -o logs/tensorqtl_mdd_amyg.txt
#$ -e logs/tensorqtl_mdd_amyg.txt
#$ -m e

#module load tensorqtl
USAGE_CUTOFF=10
NUM_GPUS=1

eval "$(/opt/conda/bin/conda shell.bash hook)"
conda activate tf

avail_gpus=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader | cut -d " " -f 1 | awk -v usage="$USAGE_CUTOFF" '$1 < usage {print NR - 1}')
echo "Available GPUs: $avail_gpus"
export CUDA_VISIBLE_DEVICES=$(echo "$avail_gpus" | head -n $NUM_GPUS | paste -sd ",")
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

./run_tensorqtl.py >& run_tensorqtl.log
