#!/bin/bash
#SBATCH --job-name=11tiss___intrinsic_subtype_2_acat_eqtl
#SBATCH --time=4:00:00
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/acat_eqtl/%x.o%j.log
### Standard Error and Output to same file

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

python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/spredixcan_expression_acat.py \
--results_dir ../../01_spredixcan_eqtl_sqtl/output/spredixcan_eqtl_mashr/ \
--study_pattern "spredixcan_igwas_gtexmashrv8_intrinsic_subtype_2*" \
--tiss_list ../input/tissue_lists/11tiss.txt \
--output ../output/acat_eqtl/11tiss__intrinsic_subtype_2_eqtl_acat_results.tsv