#!/bin/bash
#SBATCH --job-name=post_process_subtype_2
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=16GB
#SBATCH --output=logs/%x.out

cd $SLURM_SUBMIT_DIR

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/gao-lab/Julian/software/envs/metaxcan

# Replace the below arguments with subtype gwas specific details
python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_summary_imputation_postprocess.py \
-gwas_file ../output/metaxcan_gwas_parse/intrinsic_subtype_2.txt.gz \
-folder ../output/metaxcan_gwas_impute_regions/ \
-pattern intrinsic_subtype_2_* \
-parsimony 9 \
-output ../output/metaxcan_gwas_imputed/intrinsic_subtype_2.txt.gz
