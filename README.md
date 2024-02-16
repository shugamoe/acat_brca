# Overview

Each folder here has the following format:
  * `<run_order_number>_<description_of_what_is_done>`/

Each folder has the following subfolders, each with their own README.md:
  * `code/`
  * `input/`
  * `output/`

Please see the individual subfolder README.md files for further information.

## Folder summaries

`00_prep_inputs_for_metaxcan/`
  * This folder prepares raw GWAS information for MetaXcan's SPrediXcan process.
  * This includes METAL, GWAS parsing, and summary imputation to variants used in predictdb.org's GTEX V8 mashr models.

`010_ukb_bcac_meta_new_loci/`
  * This folder contains a GWAS analysis of the meta-analysis of UKB/BCAC for
    overall breast cancer to find newly identified GWAS variants.

`01_spredixcan_eqtl_sqtl/`
  * This folder contains the SPrediXcan related code/input/output.

`02_acat_eqtl_sqtl`
  * Contains the ACAT (p-value combination method) results.

`03_enloc/`
  * Finemapping enloc.

`04_cojo_ctwas/`
  * Contains conditional joint analysis (COJO) and the subsequent conditional
    TWAS (CTWAS) based on the COJO results.
  * Due to computational strain only Bonferroni significant variants from each
    study are run through the COJO/CTWAS process.

05_focus:
  * Contains Bogdan Lab's FOCUS resuls for our studies of interest.

08_expression_pred_ccv:
  * Contains code/input/output pertaining to the enrichment of credible causal variant (CCVs), made recently to address AJHG reviewer comments.

999_reporting_friendly_tables:
  * Contains precursor tables to those final tables presented in the paper.
    Convenient files for genetic distances, gencode 26 to 40 conversion, are included.
