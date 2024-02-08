#!/bin/bash
#PBS -N BCAC_Overall_BreastCancer_EUR_acat_sqtl_input
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l mem=16gb
#PBS -l nodes=1:ppn=1

#PBS -o logs/acat_prep_sqtl/${PBS_JOBNAME}.o${PBS_JOBID}.log
### Standard Error and Output to same file
#PBS -j oe

module load gcc/6.2.0
module load python/3.5.3

cd $PBS_O_WORKDIR 

export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/spredixcan_splicing_merge.py \
--spred_sqtl_pattern "../../01_spredixcan_eqtl_sqtl/output/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_BCAC_Overall_BreastCancer_EUR*.csv" \
--output ../input/acat_prep_sqtl/BCAC_Overall_BreastCancer_EUR_sqtl_acat_input.csv