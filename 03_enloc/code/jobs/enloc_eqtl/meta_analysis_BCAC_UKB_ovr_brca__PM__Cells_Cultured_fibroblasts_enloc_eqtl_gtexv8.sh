#!/bin/bash
#PBS -N meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts_enloc_eqtl_gtexv8
#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l mem=6gb
#PBS -l nodes=1:ppn=1

#PBS -o logs/enloc_eqtl/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -j oe

cd $PBS_O_WORKDIR 

if [ -f ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts.enloc.rst ]; then
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
printf "gwas_data\t../input/sliced_gwas/meta_analysis_BCAC_UKB_ovr_brca.sliced.gz\n" >> $1
printf "qtl_fm_dir\t/scratch/jmcclellan/dapg/eqtl/results/collapsed_dapg/Cells_Cultured_fibroblasts\n" >> $1
printf "out_dir\t../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts\n" >> $1
printf "trait\tmeta_analysis_BCAC_UKB_ovr_brca\n" >> $1
# bypass_enrichment_analysis 1
printf "use_openmp\t0\n" >> $1
}

[ -d ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts ] || mkdir -p ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts

create_config ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts/enloc.params

perl /gpfs/data/gao-lab/Julian/software/enloc_bin/enloc ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts/enloc.params
mv ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts/meta_analysis_BCAC_UKB_ovr_brca.enloc.rst ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts.enloc.rst
mv ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts/meta_analysis_BCAC_UKB_ovr_brca.enrich.est ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts.enrich.est
rm -rf ../output/enloc_eqtl/meta_analysis_BCAC_UKB_ovr_brca__PM__Cells_Cultured_fibroblasts

>&2 echo "Finished job."