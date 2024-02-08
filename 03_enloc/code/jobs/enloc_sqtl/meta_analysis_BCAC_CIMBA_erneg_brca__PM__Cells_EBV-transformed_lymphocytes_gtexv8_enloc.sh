#!/bin/bash
#PBS -N meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes_gtexv8_enloc
#PBS -S /bin/bash
#PBS -l walltime=48:00:00
#PBS -l mem=23gb
#PBS -l nodes=1:ppn=1

#PBS -o logs/enloc_sqtl/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -j oe

cd $PBS_O_WORKDIR 

if [ -f ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes.enloc.rst ]; then
    echo "Output present!"
    exit 0
fi


module load gcc/6.2.0
module load gsl/2.3
module load boost/1.61.0
module load bzip2/1.0.6
module load perl/5.24.0

## executable options

create_config()
{
printf "bin_dir\t/gpfs/data/gao-lab/Julian/software/enloc_bin/\n" > $1
printf "gwas_data\t../input/sliced_gwas/meta_analysis_BCAC_CIMBA_erneg_brca.sliced.gz\n" >> $1
printf "qtl_fm_dir\t/scratch/jmcclellan/dapg/sqtl/scratch/abarbeira3/v8_process/dapg/sqtl/results/enloc_ready_maf0.01_w1000000/Cells_EBV-transformed_lymphocytes\n" >> $1
printf "out_dir\t../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes\n" >> $1
printf "trait\tmeta_analysis_BCAC_CIMBA_erneg_brca\n" >> $1
# bypass_enrichment_analysis 1
printf "use_openmp\t0\n" >> $1
}

[ -d ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes ] || mkdir -p ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes

create_config ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes/enloc.params

perl /gpfs/data/gao-lab/Julian/software/enloc_bin/enloc ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes/enloc.params
mv ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes/meta_analysis_BCAC_CIMBA_erneg_brca.enloc.rst ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes.enloc.rst
mv ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes/meta_analysis_BCAC_CIMBA_erneg_brca.enrich.est ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes.enrich.est
rm -rf ../output/enloc_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca__PM__Cells_EBV-transformed_lymphocytes

>&2 echo "Finished job."