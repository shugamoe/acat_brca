#!/bin/bash
#PBS -N meta_analysis_BCAC_CIMBA_erneg_brca__PM__Stomach_spredixcan_igwas_gtexmashrv8_sqtl
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
--gwas_file ../../00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
--model_db_path ../input/gtex_v8_sqtl_dbs_mashr/mashr_Stomach.db \
--covariance ../input/gtex_v8_sqtl_dbs_mashr/mashr_Stomach.txt.gz \
--keep_non_rsid --additional_output --model_db_snp_key varID \
--throw \
--output_file ../output/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_meta_analysis_BCAC_CIMBA_erneg_brca__PM__Stomach.csv