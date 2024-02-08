#!/bin/bash
#PBS -N meta_analysis_BCAC_UKB_ovr_brca__PM__Skin_Sun_Exposed_Lower_leg_spredixcan_igwas_gtexmashrv8_eqtl
#PBS -S /bin/bash
#PBS -l walltime=0:30:00
#PBS -l mem=2gb
#PBS -l nodes=1:ppn=1

#PBS -o logs/spredixcan_eqtl_mashr/${PBS_JOBNAME}.o${PBS_JOBID}.log

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
--gwas_file ../../00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/meta_analysis_BCAC_UKB_ovr_brca.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
--model_db_path /gpfs/data/gao-lab/splicing_data_results/cp_abarbeira/projects/gtex_v8/models_v1/eqtl/mashr/mashr_Skin_Sun_Exposed_Lower_leg.db \
--covariance /gpfs/data/gao-lab/splicing_data_results/cp_abarbeira/projects/gtex_v8/models_v1/eqtl/mashr/mashr_Skin_Sun_Exposed_Lower_leg.txt.gz \
--keep_non_rsid --additional_output --model_db_snp_key varID \
--throw \
--output_file ../output/spredixcan_eqtl_mashr/spredixcan_igwas_gtexmashrv8_meta_analysis_BCAC_UKB_ovr_brca__PM__Skin_Sun_Exposed_Lower_leg.csv