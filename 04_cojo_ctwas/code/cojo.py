#!/usr/bin/env python
# -*- coding: utf-8 -*-

from parsl.app.app import bash_app, python_app
from config import id_for_memo_File
import parsl
from helpers import SetupCOJO

setup = SetupCOJO()

    # group	og_snp_list_path	imp_snp_list_path	cond_snp_list_path
@bash_app(cache=True)
def gcta_cojo(study, gene_name, snp_list_path, cond_snp_list_path, chr_num,
              out_prefix, inputs=[], outputs=[], stdout=parsl.AUTO_LOGNAME,
              stderr=parsl.AUTO_LOGNAME): 
    bfile = setup.COJO_BFILE_PATTERN.format(chr_num=chr_num)
    cojo_ma_file = inputs[0]
    output = outputs[0]
    
    extract_path = "{}{}.snplist".format(setup.COJO_IN_DIR_KEY[study], gene_name)
    pre_cojo = f"""
    IN_DIR=`dirname {extract_path}`
    mkdir -p $IN_DIR
    cat {snp_list_path} {cond_snp_list_path} | sort --unique > {extract_path}
    OUT_DIR=`dirname {out_prefix}`
    mkdir -p $OUT_DIR
    """

    bash_command = \
    f"""
    if test -f "{output}"; then
      echo "Output exists: {output}"
    else
      {pre_cojo}

      gcta-1.94.1 \
      --bfile {bfile} \
      --cojo-file {cojo_ma_file} \
      --cojo-cond {cond_snp_list_path} \
      --extract {extract_path} \
      --out {out_prefix}
    fi
    """
    return bash_command

@bash_app(cache=True)
def special_gcta_cojo(study, inputs=[], outputs=[], stdout=parsl.AUTO_LOGNAME,
              stderr=parsl.AUTO_LOGNAME): 
    
    bash_command = \
    f"""
    cp SMAD3.cond.snplist ../output/intermediate_data/clump_cond_snp_lists/{study}/ENSG00000166949.15.clumped.trimmed

    gcta-1.94.1 \
    --bfile ../input/COJO_REF_LD/jcm_hg38_panel_var_id/bfiles_merged/hg38_EUR_chr_15_gtex_v8_1000G.phase3.genotypes.final \
    --cojo-file ../output/intermediate_data/cojo_ma_files/{study}_gwas.ma \
    --cojo-cond SMAD3.cond.snplist \
    --extract SMAD3.pred.snplist \
    --out ../output/intermediate_data/cojo_output/{study}/ENSG00000166949.15
    """
    return bash_command

@python_app(cache=True)
def special_add_non_effect_allele(study, cojo_res_sep="\t", cojo_ma_sep=" ", inputs=[],
        outputs=[], stdout=parsl.AUTO_LOGNAME, stderr=parsl.AUTO_LOGNAME):
    import numpy as np
    import pandas as pd
    from glob import glob
    import os
    import re
    cojo_ma_df = pd.read_csv("../output/intermediate_data/cojo_ma_files/" + study + "_gwas.ma", sep=cojo_ma_sep)

    cojo_oput_path = "../output/intermediate_data/cojo_output/" + study + "/ENSG00000166949.15.cma.cojo"
    oput_df = pd.read_csv(cojo_oput_path, sep=cojo_res_sep)
    gene = re.split(".cma.cojo", os.path.basename(cojo_oput_path))[0] # Could be a "protein" or "sequence id" maybe soon too

    # Put in non effect allele if not present
    if 'non_effect_allele' not in oput_df:
        print(study)
        print(gene)
        oput_df['conditional_status'] = None
        # If there were any SNPs used to model the gene that overlapped with the original index SNP list, join them in
        # if necessary and zero out their effect size/standard error
        model_cond_olap_path = f"{setup.COND_SNP_LIST_DIR_KEY[study]}{gene}.snplist.model.cond.overlap"
        if os.path.exists(model_cond_olap_path):
            model_cond_olap_df = pd.read_csv(model_cond_olap_path)
            model_cond_olap_df.rename(columns={"model_snp": "SNP"}, inplace=True)
            oput_df = oput_df.merge(model_cond_olap_df, left_on='SNP', right_on='SNP', how='outer')

            # Only zero out conditional overlaps if we otherwise failed to calculate the effect
            # E.g. there is overlap, but the overlapped SNP is dropped of the condition set in clumping
            # and we were able to calculate conditional effects. So we don't throw the effect away.
            oput_df['pC'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['pC'].isnull())] = 1
            oput_df['bC'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC'].isnull())] = 0
            oput_df['conditional_status'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC_se'].isnull())] = "null_result_and_condition_overlap"
            oput_df['bC_se'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC_se'].isnull())] = 0

        # Add in any predictors in the model that were imputed. 0 them out here as well.
        # With the zscore and standard error estimation there shouldn't be any imputed SNPs remaining
        gene_all_pred_snps_path = f"{setup.STUDY_SNP_LIST_DIR_KEY[study]}{gene}.snplist.all.for.prediction"
        if os.path.exists(gene_all_pred_snps_path):
            model_all_pred_snps_df = pd.read_csv(gene_all_pred_snps_path, names=["SNP"])
            oput_df = oput_df.merge(model_all_pred_snps_df, left_on='SNP', right_on='SNP', how='outer')

        oput_df = oput_df.merge(cojo_ma_df, left_on='SNP', right_on='chromosome_position')

        # If we inserted any new variant, make sure and Chr, refA
        # (effect allele), bp, and freq are added.
        # Get chromosome
        oput_df['Chr'] = oput_df['SNP'][0].split("_")[0].split("chr")[1]
        oput_df['refA'] = oput_df['effect_allele']
        oput_df['freq'] = oput_df['frequency']

        # For any NaN conditional effects (bC, BC_se, pC), set them to 0, 0, 1
        oput_df['pC'][oput_df['pC'].isnull()] = 1
        oput_df['bC'][oput_df['bC'].isnull()] = 0
        oput_df['conditional_status'][oput_df['bC_se'].isnull()] = "freq_mismatch_or_not_in_ld_or_collinear"
        oput_df['bC_se'][oput_df['bC_se'].isnull()] = 0
        oput_df['conditional_status'][oput_df['conditional_status'].isnull()] = "cojo"

        oput_df[["throwaway_chr", "bp"]] = oput_df['SNP'].str.split("_", expand=True)
        oput_df = oput_df[['Chr', 'SNP', 'bp', 'refA', 'non_effect_allele', 'freq', 'b', 'se', 'p', 'n', 'freq_geno', 'bC', 'bC_se', 'pC', 'conditional_status']]

        new_out = f"{cojo_oput_path}"
        written = oput_df.to_csv(new_out, sep=cojo_res_sep, index=False)
        print("written")
        print(" ")
    return(written)
