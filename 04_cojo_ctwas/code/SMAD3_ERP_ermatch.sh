# SMAD3_ERP_ermatch.sh
export PATH=/gpfs/data/huo-lab/jmcclellan/software/gcta-1.94.1-linux-kernel-3-x86_64:$PATH;

gcta-1.94.1 \
--bfile ../input/COJO_REF_LD/jcm_hg38_panel_var_id/bfiles_merged/hg38_EUR_chr_15_gtex_v8_1000G.phase3.genotypes.final \
--cojo-file ../output/intermediate_data/cojo_ma_files/BCAC_ERPOS_BreastCancer_EUR_index_ermatch_gwas.ma \
--cojo-cond SMAD3.cond.snplist \
--extract SMAD3.pred.snplist \
--out ../output/intermediate_data/cojo_output/BCAC_ERPOS_BreastCancer_EUR_index_ermatch/ENSG00000166949.15
