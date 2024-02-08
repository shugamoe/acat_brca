#!/bin/bash
#SBATCH --job-name=11tiss___meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_results
#SBATCH --time=4:00:00
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/acat_sqtl/%x.o%j.log

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

python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/SMulTiXcanByFeature.py \
--gwas_file ../../00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/meta_analysis_BCAC_UKB_ovr_brca.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore  \
--grouping ../input/combine_phenotype_groups.txt.gz GTEx_sQTL \
--model_db_path ../input/combine_sqtl.db \
--keep_non_rsid --model_db_snp_key varID \
--covariance ../input/combine_cross_intron_covar.txt.gz \
--associations ../input/acat_prep_sqtl/meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_input.csv SPrediXcan \
--cutoff_condition_number 30 \
--verbosity 10 \
--acat \
--tiss_list ../input/tissue_lists/11tiss.txt \
--output ../output/acat_sqtl/11tiss__meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_results.tsv