#!/bin/bash
#SBATCH --job-name=dz_jl_tables
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=32GB
#SBATCH --output=logs/%x.out
#SBATCH --partition=tier2q

cd $SLURM_SUBMIT_DIR

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/huo-lab/jmcclellan/software/envs/parsl

Rscript dz_jl_tables_for_jcm.R
