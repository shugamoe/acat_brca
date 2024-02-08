#!/bin/bash
#PBS -N slice_BCAC_Overall_BreastCancer_EUR
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l mem=4gb
#PBS -l nodes=1:ppn=1
#PBS -o logs/slice_gwas/${PBS_JOBNAME}.o${PBS_JOBID}.log
### Standard Error and Output to same file
#PBS -j oe
########################################################################################################################
#CRI submission dandruff

module load gcc/6.2.0
module load python/3.5.3

export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/slice_gwas_by_region.py \
-region_file ../../00_prep_inputs_for_metaxcan/input/eur_ld.bed.gz \
-gwas_file ../../00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/BCAC_Overall_BreastCancer_EUR.txt.gz \
-parsimony 8 \
-output ../input/sliced_gwas/BCAC_Overall_BreastCancer_EUR.sliced.gz