#!/bin/bash
#SBATCH --job-name=intrinsic_subtype_4_acat_sqtl_input
#SBATCH --time=4:00:00
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/acat_prep_sqtl/%x.o%j.log

# gardner
# module load gcc/6.2.0
# module load python/3.5.3

# randi
module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/gao-lab/Julian/software/envs/metaxcan
export PATH=/gpfs/data/gao-lab/Julian/software/envs/metaxcan/bin:$PATH


cd $SLURM_SUBMIT_DIR 

export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/spredixcan_splicing_merge.py \
--spred_sqtl_pattern "../../01_spredixcan_eqtl_sqtl/output/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_intrinsic_subtype_4*.csv" \
--output ../input/acat_prep_sqtl/intrinsic_subtype_4_sqtl_acat_input.csv