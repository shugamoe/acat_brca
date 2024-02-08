#!/bin/bash
#SBATCH --job-name=parse_Luminal_B_HER2Neg
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=32GB
#SBATCH --output=logs/%x.out

cd $SLURM_SUBMIT_DIR

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/gao-lab/Julian/software/envs/metaxcan


# Replace the below arguments with subtype gwas specific details
python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_parsing.py \
-gwas_file ../input/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt \
-output_column_map var_name variant_id -output_column_map Baseline.Meta non_effect_allele -output_column_map Effect.Meta effect_allele \
-output_column_map Luminal_B_HER2Neg_log_or_meta effect_size -output_column_map Luminal_B_HER2Neg_se_meta standard_error -output_column_map Freq1 frequency \
-split_column variant_id '_' chromosome position nea ea --chromosome_format \
--insert_value sample_size 319209 --insert_value n_cases 73067 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-liftover ../input/hg19ToHg38.over.chain.gz \
-snp_reference_metadata ../input/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
-output ../output/metaxcan_gwas_parse/intrinsic_subtype_Luminal_B_HER2Neg.txt.gz
