#!/bin/bash
#SBATCH --job-name=intrinsic_subtype_2__PM__Minor_Salivary_Gland_spredixcan_igwas_gtexmashrv8_sqtl
#SBATCH --time=0:30:00
#SBATCH --mem=2gb
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/spredixcan_sqtl_mashr/%x.o%j.log

# gardner
# module load gcc/6.2.0
# module load python/3.5.3

# randi
module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/gao-lab/Julian/software/envs/metaxcan
export PATH=/gpfs/data/gao-lab/Julian/software/envs/metaxcan/bin:$PATH


export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd $SLURM_SUBMIT_DIR 

python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/SPrediXcan.py \
--gwas_file ../../00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/intrinsic_subtype_2.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
--model_db_path ../input/gtex_v8_sqtl_dbs_mashr/mashr_Minor_Salivary_Gland.db \
--covariance ../input/gtex_v8_sqtl_dbs_mashr/mashr_Minor_Salivary_Gland.txt.gz \
--keep_non_rsid --additional_output --model_db_snp_key varID \
--throw \
--output_file ../output/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_intrinsic_subtype_2__PM__Minor_Salivary_Gland.csv