#!/bin/bash
#PBS -N meta_analysis_BCAC_UKB_ovr_brca__parse_e
#PBS -S /bin/bash
#PBS -l walltime=6:00:00
#PBS -l mem=4gb
#PBS -l nodes=1:ppn=1

#PBS -o logs/parse_enloc_result_sqtl/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -j oe

module load gcc/6.2.0
module load python/3.5.3

cd $PBS_O_WORKDIR 

python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/post_process_and_merge.py \
-input_folder ../output/parse_enloc_result_sqtl \
-input_pattern "(.*)__PM__(.*).enloc.rst.gz" \
-name_subfield trait 1 -name_subfield tissue 2 \
--input_filter meta_analysis_BCAC_UKB_ovr_brca__PM__.*enloc.rst.gz \
-trait_map public_gwas_filename_map.yaml \
-output ../output/final_enloc_result_sqtl/meta_analysis_BCAC_UKB_ovr_brca__parse_e.txt