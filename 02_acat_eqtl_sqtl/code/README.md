010_acat_eqtl.jinja
010_acat_eqtl.sh
010_acat_eqtl.yaml
  * `source 010*.sh` runs expression based ACAT for all studies, and for all tissue lists.
  * (Tissue lists are either 11 tissues or "all" GTEx V8 tissues.)

020_acat_prep_sqtl.jinja
020_acat_prep_sqtl.sh
020_acat_prep_sqtl.yaml
  * `source 020*.sh` combines previous spredixcan sqtl output into a single file per study.
  * Output is stored in `input/acat_prep_sqtl`.

021_merge_sqtl_db.py
  * Only has to be run once, combines all sqtl model tissue db files into a
  * single db to support multi-tissue ACAT, preserves tissue of origin in the
  * model.

022_merge_phenotype_groups.py
  * Only has to be run once, combines all intron to ENSG gene mappings at a per
  * tissue level to a single file.

025_acat_tiss_breakdown_sqtl.jinja
025_acat_tiss_breakdown_sqtl.yaml
025_run_acat_tiss_breakdown_sqtl.sh
  * `source 025*.sh` runs splicing based ACAT for all studies but formats the
    results so that within-tissue ACAT results for a gene's introns are easily
    accessible, as well as individual intron p-values.

029_acat_sqtl.jinja
029_acat_sqtl.yaml
029_run_acat_sqtl.sh
  * `source 029*.sh` runs splicing based ACAT for all studies, and for all tissue lists.
  * (Tissue lists are either 11 tissues or "all" GTEx V8 tissues.)

jobs
logs
