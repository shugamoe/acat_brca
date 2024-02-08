#!/bin/bash
#PBS -N BCAC_Overall_BreastCancer_EUR__PM__Cells_Cultured_fibroblasts_spredixcan_igwas_gtexmashrv8_sqtl
#PBS -S /bin/bash
#PBS -l walltime=0:30:00
#PBS -l mem=2gb
#PBS -l nodes=1:ppn=1

#PBS -o logs/spredixcan_sqtl_mashr/${PBS_JOBNAME}.o${PBS_JOBID}.log

### Standard Error and Output to same file
#PBS -j oe

module load gcc/6.2.0
module load python/3.5.3

export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR 

python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/SPrediXcan.py \
--gwas_file ../../00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/BCAC_Overall_BreastCancer_EUR.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
--model_db_path /gpfs/data/gao-lab/splicing_data_results/cp_abarbeira/projects/gtex_v8/models_v1/sqtl/mashr/mashr_Cells_Cultured_fibroblasts.db \
--covariance /gpfs/data/gao-lab/splicing_data_results/cp_abarbeira/projects/gtex_v8/models_v1/sqtl/mashr/mashr_Cells_Cultured_fibroblasts.txt.gz \
--keep_non_rsid --additional_output --model_db_snp_key varID \
--throw \
--output_file ../output/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_BCAC_Overall_BreastCancer_EUR__PM__Cells_Cultured_fibroblasts.csv