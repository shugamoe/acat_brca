#!/bin/bash
#PBS -N BCAC_Overall_BreastCancer_EUR_chr8_sb1_by_region
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=1
#PBS -o logs/summary_imputation/${PBS_JOBNAME}.o${PBS_JOBID}.log
### Standard Error and Output to same file
#PBS -j oe
########################################################################################################################
#CRI submission dandruff

module load gcc/6.2.0
module load python/3.5.3

export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_summary_imputation.py \
-by_region_file ../input/eur_ld.bed.gz \
-gwas_file ../output/metaxcan_gwas_parse/BCAC_Overall_BreastCancer_EUR.txt.gz \
-parquet_genotype ../input/gtex_v8_parquet_eur_maf0.01_biallelic/gtex_v8_eur_itm.chr8.variants.parquet \
-parquet_genotype_metadata ../input/gtex_v8_parquet_eur_maf0.01_biallelic/gtex_v8_eur_itm.variants_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 8 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 1 \
--standardise_dosages \
-output ../output/metaxcan_gwas_impute_regions/BCAC_Overall_BreastCancer_EUR_chr8_sb1_reg0.1_ff0.01_by_region.txt.gz