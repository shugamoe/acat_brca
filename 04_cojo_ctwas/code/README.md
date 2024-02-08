## COJO

config.py
  * Denotes the university computing resources for use in the COJO part of the COJO/CTAS process.
  * Expects TORQUE, can be switched to slurm or other resource managers (or local resources only).
  * For more information, see the documentation for the [python library: parsl.](https://parsl.readthedocs.io/en/stable/)

cojo.py
  * A separate file that has only one function for use with the parsl library, the final gcta calculation function.

helpers.py
  * Big setup file that connects GWAS study information, e.g. sample size, file locations, genes to run etc.

prepare_data.py
  * File that contains most of the non gcta_cojo functions used in the COJO process.

run.py
  * `python3 run.py` Runs the COJO part of COJO/CTWAS, for all studies, expression and splicing.
  * Clumping R2 in use is `0.65`
  * We look 2MB left and right of a gene stat/end for previously reported GWAS variants.
  * We include ~10 or fewer variants we discovered in the BCAC/UKB METAL GWAS in the above GWAS variants.


## CTWAS

ctwas:
  * This folder contains conditional TWAS section of the COJO/CTWAS analysis. It relies on output from COJO. 
  
  * In the ctwas folder, you will find folders of the form:
    * "{study_name}_{'expression' if expression else if splicing ''}"
    00_prep_inputs_for_metaxcan/code

  * All folders, regardless of whether they are expression or splicing, contains the following files:
  |-- 00_gwas_parsing.jinja
  |-- 00_gwas_parsing.sh
  |-- 00_gwas_parsing.yaml
    * These files parse the COJO results, much like `00_prep_inputs_for_metaxcan/code/*parse*`
    * A key difference is that instead of having one result file per study, we have one result per gene of interest.


  |-- 031_conditional_covar_by_tiss.jinja
  |-- 031_conditional_covar_by_tiss.sh
  |-- 031_conditional_covar_by_tiss.yaml
  |-- 031_find_gene_snps_cond_and_pred.sh -> ../utils/031_find_gene_snps_cond_and_pred.sh
    * These files look back to the relevant COJO intermediate data for the study,
    * expression or splicing, and calculates new individual tissue conditional covariance matrixes
    * using the same conditioning and modelling variants used in the prior COJO process.

  |-- 03_metaxcan.jinja1_spredixcan_eqtl_sqtl/code
  |-- 03_metaxcan.sh
    * These files  are similar to those found in `1_spredixcan_eqtl_sqtl/code`, but final results files are on a per gene/intron per tissue basis, not a study per tissue basis.

  |-- 03_metaxcan_run_staggered_eqtl.py -> ../utils/03_metaxcan_run_staggered_eqtl.py
  |-- 03_metaxcan_run_staggered_sqtl.py -> ../utils/03_metaxcan_run_staggered_sqtl.py
  |-- 03_metaxcan_eqtl_mashr.yaml
  |-- 03_metaxcan_sqtl_mashr.yaml
    * Depending on whether a folder is for making ctwas eqtl or sqtl results, 2 relevant files from above will be found.
### CTWAS Splicing
  * If a folder is a "splicing" folder ("{study_name}")
  |-- 05_multi_tissue_intronxcan.jinja
  |-- 05_multi_tissue_intronxcan.sh
  |-- 05_multi_tissue_intronxcan.yaml
    * These files above are similar to those in `02_acat_eqtl_sqtl/code/029*`,
    * but they output results by gene/tissue list (either all or 11 tissues),
    * rather than by study/tissue list.

  |-- 051_combine_compare_mtiss_ixcan.py -> ../utils/051_combine_compare_mtiss_ixcan.py
    * This combines the above files into 1 file per study per tissue list.
    * Comparison is no longer done in this file, as it introduces some weird
    * circularity into the analysis that we can avoid.
  
### CTWAS Expression
  * If a folder is a "splicing" folder ("{study_name}_expression")
  |-- 06_1step_acat_and_combine.py -> ../utils/06_1step_acat_and_combine.py
    * This file calculates acat the same way as in
    * `02_acat_eqtl_sqtl/code/010_acat_eqtl*`, but is leaner code that doesn't
    * rely on MetaXCan code.