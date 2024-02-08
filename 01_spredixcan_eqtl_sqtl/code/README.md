For more information on running S-PrediXcan, see the original [github wiki](https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS#s-predixcan)

010_metaxcan.jinja
010_metaxcan_imputed_gtexv8_eqtl_mashr.yaml
010_run_metaxcan_imputed_gtexv8_eqtl_mashr.sh
  * `source 010*.sh` Runs GTEX V8 eqtl mashr models with all studies.

020_metaxcan.jinja
020_metaxcan_imputed_gtexv8_sqtl_mashr.yaml
020_run_metaxcan_imputed_gtexv8_sqtl_mashr.sh
  * `source 020*.sh` Runs GTEX V8 sqtl mashr models with all studies.

jobs:
spredixcan_eqtl_mashr  spredixcan_sqtl_mashr

logs
spredixcan_eqtl_mashr  spredixcan_sqtl_mashr