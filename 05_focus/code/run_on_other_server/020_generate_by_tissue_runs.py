#!/software/python-anaconda-2020.02-el7-x86_64/bin/python
TISSUES = [
'adipose_subcutaneous',
'adipose_visceral_omentum',
'breast_mammary_tissue',
'cells_cultured_fibroblasts',
'cells_ebv-transformed_lymphocytes',
'liver',
'ovary',
'spleen',
'uterus',
'vagina',
'whole_blood'
]

PATTERN="""#!/bin/bash
#SBATCH --job-name=f_{eqtl_sqtl}
#SBATCH --array=1-22%4
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH	--nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --tasks-per-node=1
#SBATCH --mem=16GB
#SBATCH --output=logs/focus_{eqtl_sqtl}_{cur_tissue}.out

module load python
module load pigz
source activate /project2/guiming/Julian/software/envs/focus

# export PATH=`python -m site --user-base`/bin/:$PATH

i=$SLURM_ARRAY_TASK_ID
cd $SLURM_SUBMIT_DIR
STUDY_NAME=$(basename `readlink -f ../`)

focus finemap "../input/munged_${{STUDY_NAME}}.sumstats.gz" \
  "../input/by_chrom/gtex_v8_eur_filtered_has_missing_snps_chr${{i}}" \
  "../input/combined_dbs/{cur_tissue}_{eqtl_sqtl}.db" \
  --verbose \
  --tissue {cur_tissue} \
  --locations ../input/regions/chr${{i}}.txt \
  --p-threshold 1e-3 \
  --chr ${{i}} \
  --out "../output/by_tissue/{eqtl_sqtl}_focus_11tiss_${{STUDY_NAME}}_Chr${{i}}__PM__{cur_tissue}"

"""

if __name__ == '__main__':
    for cur_tissue in TISSUES:
        for cur_eqtl_sqtl in ['eqtl', 'sqtl']:
            cur_pattern = PATTERN.format(cur_tissue=cur_tissue,
              eqtl_sqtl=cur_eqtl_sqtl) 
            with open('by_tissue_runs/08_focus_{}_{}.sbatch'.format(cur_eqtl_sqtl, cur_tissue), 'w') as f:
                f.write(cur_pattern)
