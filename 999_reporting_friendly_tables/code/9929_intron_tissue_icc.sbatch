#!/bin/bash
#SBATCH --job-name=11tiss_intron_icc
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=64GB
#SBATCH --output=logs/%x.out

#######################################################################################################################
#CRI submission dandruff

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/huo-lab/jmcclellan/software/envs/parsl

cd $SLURM_SUBMIT_DIR

Rscript 9929_intron_tissue_icc.R
