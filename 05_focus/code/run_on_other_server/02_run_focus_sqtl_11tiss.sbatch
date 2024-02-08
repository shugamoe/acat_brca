#!/bin/bash
#SBATCH --job-name=f_sqtl
#SBATCH --array=1-22%4
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH	--nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --tasks-per-node=1
#SBATCH --mem=16GB
#SBATCH --output=logs/focus_sqtl%a.out

module load python
module load pigz
source activate /project2/guiming/Julian/software/envs/focus

# export PATH=`python -m site --user-base`/bin/:$PATH

i="${SLURM_ARRAY_TASK_ID}"
cd $SLURM_SUBMIT_DIR
STUDY_NAME=$(basename `readlink -f ../`)

focus finemap "../input/munged_${STUDY_NAME}.sumstats.gz" \
  "../input/by_chrom/gtex_v8_eur_filtered_has_missing_snps_chr${i}" \
  "../input/combined_dbs/11tiss_sqtl.db" \
  --verbose \
  --tissue breast_mammary_tissue \
  --locations ../input/regions/chr${i}.txt \
  --p-threshold 1e-3 \
  --chr ${i} \
  --out "../output/sqtl_focus_11tiss_${STUDY_NAME}_Chr${i}"
