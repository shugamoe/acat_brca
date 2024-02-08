#!/bin/bash
#SBATCH --job-name=gene_rollup_ctwas_eqtl
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=32GB
#SBATCH --partition=tier2q
#SBATCH --output=logs/%x.out

cd $SLURM_SUBMIT_DIR

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/gao-lab/Julian/software/envs/parsl

Rscript gene_rollup_ctwas_eqtl.R
