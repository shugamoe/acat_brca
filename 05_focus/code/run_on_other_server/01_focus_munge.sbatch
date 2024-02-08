#!/bin/bash
#SBATCH --job-name=f_munge
#SBATCH --time=5:00:00
#SBATCH --partition=broadwl
#SBATCH	--nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=16GB
#SBATCH --output=logs/focus_munge.out

module load python
module load pigz
# conda activate /project2/guiming/Julian/software/envs/focus
source activate /project2/guiming/Julian/software/envs/focus

cd $SLURM_SUBMIT_DIR
STUDY_NAME=$(basename `readlink -f ../`)

# Case/control numbers double checked by Guimin, used to calculate effective sample sizes for COJO/CTWAS
focus munge "../input/${STUDY_NAME}.txt.gz" \
  --N 424601 \
  --N-cas 133511 \
  --output "../input/munged_${STUDY_NAME}"
