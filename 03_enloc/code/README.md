010_run_slice_gwas.sh
010_slice_gwas.jinja
010_slice_gwas.yaml
  * Uses hakyimlab's `summary-gwas-imputation`
  * [repo](https://github.com/shugamoe/summary-gwas-imputation) to slice the
  * GWAS for each study by region.

020_enloc.jinja
020_enloc_eqtl.yaml
020_enloc_eqtl_tissue_spec.yaml
020_run_enloc_eqtl.sh
  * `source 020*.sh` runs ENLOC for expression. Relies on pre-calculated dapg
  * results that are very large. Available on request.

021_parse_enloc_result_eqtl.jinja
021_parse_enloc_result_eqtl.yaml
021_run_parse_enloc_result_eqtl.sh
  * `source 021*.sh` runs some initial parsing for enloc eqtl results.

022_postprocess_enloc_eqtl.yaml
022_run_postprocess_enloc_eqtl.sh
  * `source 021*.sh` runs final post-processing for enloc eqtl results.

030_enloc.jinja
030_enloc_sqtl.yaml
030_enloc_sqtl_tissue_spec.yaml
030_run_enloc_sqtl.sh
  * `source 030*.sh` runs ENLOC for splicing. Relies on pre-calculated dapg
  * results that are very large. Available on request.

031_parse_enloc_result_sqtl.jinja
031_parse_enloc_result_sqtl.yaml
031_run_parse_enloc_result_sqtl.sh
  * `source 031*.sh` runs some initial parsing for enloc sqtl results.

032_postprocess_enloc_sqtl.yaml
032_run_postprocess_enloc_sqtl.sh

postprocess_and_merge.jinja
  * jinja templates file used by the 020*/030* [badger](https://github.com/hakyimlab/badger) files above.
public_gwas_filename_map.yaml
  * For use in badger, used to denote our studies of interest. 

jobs
logs
