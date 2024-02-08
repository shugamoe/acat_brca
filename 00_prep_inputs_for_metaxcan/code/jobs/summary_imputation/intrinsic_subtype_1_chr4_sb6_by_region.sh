#!/bin/bash
#SBATCH --job-name=intrinsic_subtype_1_chr4_sb6_by_region
#SBATCH --time=4:00:00
#SBATCH --mem=6gb
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/summary_imputation/%x.o%j.log
########################################################################################################################
# gardner
# module load gcc/6.2.0
# module load python/3.5.3

# randi
module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/gao-lab/Julian/software/envs/metaxcan
export PATH=/gpfs/data/gao-lab/Julian/software/envs/metaxcan/bin:$PATH

echo `which python3`
echo `module list`

export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd $SLURM_SUBMIT_DIR

python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_summary_imputation.py \
-by_region_file ../input/eur_ld.bed.gz \
-gwas_file ../output/metaxcan_gwas_parse/intrinsic_subtype_1.txt.gz \
-parquet_genotype ../input/gtex_v8_parquet_eur_maf0.01_biallelic/gtex_v8_eur_itm.chr4.variants.parquet \
-parquet_genotype_metadata ../input/gtex_v8_parquet_eur_maf0.01_biallelic/gtex_v8_eur_itm.variants_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 4 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 6 \
--standardise_dosages \
-output ../output/metaxcan_gwas_impute_regions/intrinsic_subtype_1_chr4_sb6_reg0.1_ff0.01_by_region.txt.gz