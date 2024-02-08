#!/bin/bash
#PBS -N BCAC_ERPOS_BreastCancer_EUR__parse_e
#PBS -S /bin/bash
#PBS -l walltime=6:00:00
#PBS -l mem=4gb
#PBS -l nodes=1:ppn=1

#PBS -o logs/parse_enloc_result_eqtl/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -j oe

module load gcc/6.2.0
module load python/3.5.3

cd $PBS_O_WORKDIR 

python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/post_process_and_merge.py \
-input_folder ../output/parse_enloc_result_eqtl \
-input_pattern "(.*)__PM__(.*).enloc.rst.gz" \
-name_subfield trait 1 -name_subfield tissue 2 \
--input_filter BCAC_ERPOS_BreastCancer_EUR__PM__.*enloc.rst.gz \
-trait_map public_gwas_filename_map.yaml \
-output ../output/final_enloc_result_eqtl/BCAC_ERPOS_BreastCancer_EUR__parse_e.txt