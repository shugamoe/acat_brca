For more information on running FOCUS, please see the original Github repo's [wiki.](https://github.com/bogdanlab/focus/wiki)

Our paper uses a [slightly modified version of FOCUS](https://github.com/shugamoe/focus) with changes (some from myself, others from [Alvaro Barbeira](https://github.com/heroico)) to accomodate its use with
GTEx V8 models from predictdb.org.


## Files 

00_format_gwas.pbs
00_format_gwas.py
  * These files format the relevant GWAS files from `../../00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/` into a format more suitable for FOCUS's built in "munging" (GWAS formatting) function.

run_on_other_server/
  * I had a never ending series of issues trying to install FOCUS on our [main * computing * resource.](https://wiki.uchicago.edu/display/public/CRI/Gardner+Quick+Start+Guide).
  * So I ran most FOCUS-related code on an [alternate computing
  * resource](https://rcc.uchicago.edu/).

|-- 00_patch_dbs.py
  * This file was run to replace the rsID column with the varID column in the
  * predictdb.org databases, essentially a hack to accomodate the fact that we
  * had our LD data .bim files in ~ <chr_num>_<pos>_<a1>_<a2> sort of format.

|-- 01_focus_munge.sbatch
  * Multiple versions of this file existed, depending on the case/control count
  * of the study. This is FOCUS's built in GWAS processing step.

|-- 020_generate_by_tissue_runs.py
  * Generates runs of FOCUS, 1 run per chromosome per (11) tissues.

|-- 02_run_focus_eqtl_11tiss.sbatch
|-- 02_run_focus_sqtl_11tiss.sbatch
  * Submits the eqtl/sqtl based FOCUS runs in a way that doesn't blitz SLURM.

|-- 09_coalesce_by_tissue_output.py
  * Coallesce the individual outputs by chromosome/tissue into
  * a single file per study for convenience.

|-- 2023_focus_eqtl_11tiss_one_by_one.db.sh
|-- 2023_focus_sqtl_11tiss_one_by_one.db.sh
  * These files create the pre-requisite `.db` files needed to run FOCUS. A
  * single `.db` file is created for each (11) tissue for both expression and
  * splicing.

logs
