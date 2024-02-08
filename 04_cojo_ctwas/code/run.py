#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import re
import pandas as pd

import parsl
from parsl.dataflow.futures import AppFuture
from config import config, id_for_memo_File
from parsl.data_provider.files import File
from parsl.app.app import bash_app 
from parsl.app.errors import MissingOutputs, BashExitFailure
from functools import partial
from math import ceil
from datetime import datetime

from helpers import SetupCOJO
from prepare_data import make_study_snp_lists, \
    make_cojo_ma_file, \
    make_cond_snp_lists, \
    make_cojo_run_data, \
    make_clump_cond_snp_list, \
    make_clump_og_proxy_cond_snp_list, \
    add_non_effect_allele, \
    make_alt_snp_lists

from condTWAS import parse_gene_gwas, \
    make_pred_and_cond_snps_gene_gwas, \
    make_want_genes, \
    calc_tiss_cond_covar, \
    spredixcan, \
    combine_eqtl_and_calc_acat, \
    calc_mtiss_acat, \
    combine_sqtl_acat

from cojo import gcta_cojo, special_gcta_cojo, special_add_non_effect_allele

setup = SetupCOJO()
config.retries = 0
parsl.load(config)
parsl.set_stream_logger() # log lines to screen


def run_cojo(study, index_snp_window=5000000, na_snp_window=6000,
        og_check_window=6000, expression=False, clump_r2_imputed=.7,
        collin_prevent_mode="plink", ermatch=None,
        
        run_ctwas_after=True,

        # Passed to run_ctwas
        tiss_file="../input/tissue_lists/11tiss.txt",
        has_breast=True
        ):
    genes_to_run, colname_for_genes = setup.GENES_TO_RUN[study]

    cojo_ma_file_res = make_cojo_ma_file(study, outputs=[File(os.path.join(setup.COJO_MA_DIR, f"{study}_gwas.ma"))])

    study_snp_lists_res = make_study_snp_lists(study, gene_whitelist_colname=colname_for_genes,
                                                  inputs=[cojo_ma_file_res.outputs[0], genes_to_run], outputs=[File(os.path.join(setup.STUDY_SNP_LIST_PARENT, f"{study}_tracker.tsv"))],
                               expression=expression)

    cond_snp_lists_res = make_cond_snp_lists(study, bp_range=index_snp_window,
            inputs=[study_snp_lists_res.outputs[0], cojo_ma_file_res.outputs[0]],
                                                   outputs=[File(os.path.join(setup.COND_SNP_LIST_PARENT,
                                                                              f"{study}_tracker.tsv"))], ermatch=ermatch)

    cond_snp_list  = pd.read_csv(cond_snp_lists_res.outputs[0].result(), sep="\t")

    clump_cond_snp_lists = []
    for gene in cond_snp_list.itertuples():
        if gene.status == "has_original_index_snps":
            clump_cond_snp_lists.append(
                    make_clump_cond_snp_list(study, chr_num=gene.chr_num,
                        plink_output =
                        os.path.join(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study],
                            f"{gene.group}.clumped"),
                        clump_r2 = clump_r2_imputed, # Getting collinearity probem with default .85 that worked before
                        cond_snp_list=File(gene.cond_snp_list_path),
                        cojo_ma_path=cojo_ma_file_res.outputs[0],
                        collin_prevent_mode=collin_prevent_mode,
                        outputs=[File(os.path.join(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study],
                            f"{gene.group}.clumped.trimmed"))])
            )
        elif gene.status == "has_original_index_snps_with_proxies":
            clump_cond_snp_lists.append(
                    make_clump_cond_snp_list(study, chr_num=gene.chr_num,
                        plink_output =
                        os.path.join(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study],
                            f"{gene.group}.clumped"),
                        clump_r2 = clump_r2_imputed, # Getting collinearity probem with default .85 that worked before
                        cond_snp_list=File(gene.cond_snp_list_path),
                        cojo_ma_path=cojo_ma_file_res.outputs[0],
                        collin_prevent_mode=collin_prevent_mode,
                        outputs=[File(os.path.join(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study],
                            f"{gene.group}.clumped.trimmed"))])
            )
        # elif gene.status.startswith("has_index_snps_imputed"):
        #     clump_cond_snp_lists.append(
        #             make_clump_og_proxy_cond_snp_list(study, chr_num=gene.chr_num,
        #                 plink_output =
        #                 os.path.join(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study],
        #                     f"{gene.group}.clumped.best"),
        #                 cond_snp_list=File(gene.cond_snp_list_path),
        #                 nearby_og_snp_list=File(gene.nearby_og_snps_path),
        #                 outputs=[File(os.path.join(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study],
        #                     f"{gene.group}.clumped.trimmed"))])
        #     )
        else:
            clump_cond_snp_lists.append(None)

    # Catch where clumping failed by enumerating through clump_cond_snp_lists?
    cojo_run_dat_res = make_cojo_run_data(study, study_list=study_snp_lists_res.outputs[0],
                                          cond_genes_list=cond_snp_lists_res.outputs[0],
                                          inputs = [elmnt.outputs[0] if type(elmnt)
                                              is AppFuture else None for
                                              elmnt in clump_cond_snp_lists],
                                          outputs=[File(os.path.join(setup.COJO_IN_PARENT, f"{study}_blacklist.tsv")),
                                                   File(os.path.join(setup.COJO_IN_PARENT, f"{study}_whitelist.tsv"))])
    
    # Kick off cojo runs
    gene_whitelist  = pd.read_csv(cojo_run_dat_res.outputs[1].result(), sep="\t")

    # Add in SMAD3 gene for James since he might want it and run COJO for it
    if study in ["BCAC_ERPOS_BreastCancer_EUR_index_ermatch",
            "BCAC_ERPOS_BreastCancer_EUR_index_erswap", 
            "BCAC_ERPOS_BreastCancer_EUR",
            "meta_analysis_BCAC_UKB_ovr_brca"
            ]:

        # meta_analysis_BCAC_UKB_ovr_brca has the file but still doesn't run, overwrite
        # Check if this is still needed before running it blindly, only need to
        # run if we didn't have anything to condition this gene on
        # if not os.path.exists(f"../output/intermediate_data/clump_cond_snp_lists/{study}/ENSG00000166949.15.clumped.trimmed"):
        special_res = special_gcta_cojo(study, outputs=[File(f"../output/intermediate_data/cojo_output/{study}/ENSG00000166949.15.cma.cojo")])
        special_add_non_effect_allele(study=study, inputs=[special_res])
    
    # COJO Runs
    cojo_all_res = []
    for gene in gene_whitelist.itertuples():
        out_prefix = "{}{}".format(setup.COJO_OUT_DIR_KEY[study], gene.group)
        cojo_all_res.append(gcta_cojo(study, gene.group, gene.snp_list_path,
                  gene.clump_cond_snp_list_path, gene.chr_num, out_prefix,
                  inputs=[cojo_ma_file_res.outputs[0]],
                  outputs=[File(f"{out_prefix}.cma.cojo")])
                  )


    ne_res = []
    for cojo_res in cojo_all_res: 
        ne_res.append(
        add_non_effect_allele(study, cojo_ma_path=cojo_ma_file_res.outputs[0], inputs=[cojo_res.outputs[0]]))

    # might not need to wait for current task below since add_non_effect allele now returns an app future.
    # parsl.wait_for_current_tasks()

    if run_ctwas_after is True:
        # finished_ne_res = [res.result() for res in ne_res]
        run_ctwas(study, tiss_file=tiss_file, has_breast=has_breast, cojo_final_res=ne_res)

    # No longer needed
    # In some cases we might have {gene}.cma.allna files, indicating that all
    # the genes we wanted conditional effects for could not be calculated
    # (usually when the # of genes we want is small), so we go in and "fix" 
    # these genes by creating a new list of SNPs, taking the nearest
    # <na_snp_window> original SNPs on each side of the NA SNPs and running 
    # COJO with those instead.
    # make_alt_snp_lists(study, na_snp_window, whitelist_path=cojo_run_dat_res.outputs[1].result().filepath,
    #         cojo_ma_path=cojo_ma_file_res.outputs[0].result().filepath)
    # parsl.wait_for_current_tasks()
    # gene_whitelist  = pd.read_csv(cojo_run_dat_res.outputs[1].result(), sep="\t")
    # for gene in gene_whitelist.itertuples():
    #     if gene.status == "has_original_index_snps_na_fix":
    #         out_prefix = "{}{}".format(setup.COJO_OUT_DIR_KEY[study], gene.group)
    #         gcta_cojo(study, gene.group, gene.na_fix_snp_list_path,
    #                   gene.clump_cond_snp_list_path, gene.chr_num, out_prefix,
    #                   inputs=[cojo_ma_file_res.outputs[0]],
    #                   outputs=[File(f"{out_prefix}.cma.cojo")])
    # parsl.wait_for_current_tasks()
    # add_non_effect_allele(study, cojo_ma_path=cojo_ma_file_res.outputs[0])
    return

def run_ctwas(study, tiss_file="../input/tissue_lists/11tiss.txt",
        has_breast=True, cojo_final_res=[]):
    """
    """
    from glob import glob
    from pathlib import Path
    # Gather cojo_output files, extract gene names, parse each gene output
    cojo_out_files = glob(f"../output/intermediate_data/cojo_output/{study}/*.cma.cojo")

    cojo_out_gene_names = [re.match("(ENSG\d{1,}\.\d{1,})\.cma\.cojo", thing)[1] for thing in map(os.path.basename, cojo_out_files)]


    gene_parsed_res = []
    for gene_file, gene_name in zip(cojo_out_files, cojo_out_gene_names):
        out_file = File(f"../output/intermediate_data/condTWAS_parse/{study}/{gene_name}.txt.gz")
        gene_parsed_res.append(parse_gene_gwas(gene_file, gene_name, study=study, inputs=cojo_final_res, outputs=[out_file]))
    parsl.wait_for_current_tasks()

    # For each gene file, prepare the prediction snps and condition snps to calculate conditioned covariances
    pred_cond_snps_res = []
    for gene_parse_res, gene_name in zip(gene_parsed_res, cojo_out_gene_names):
        pred_cond_snps_res.append(make_pred_and_cond_snps_gene_gwas(gene_name=gene_name, study=study, inputs=[gene_parse_res.outputs[0]],
                outputs=[File(f"../output/intermediate_data/condTWAS_cond_covar/{study}/{gene_name}.snps.pred.txt"),
                         File(f"../output/intermediate_data/condTWAS_cond_covar/{study}/{gene_name}.snps.condition.txt")])
                )
    parsl.wait_for_current_tasks()

    # Make the list of want genes
    want_genes = make_want_genes(study, inputs=[res.outputs[0] for res in gene_parsed_res],
            outputs=[File(f"../output/intermediate_data/condTWAS_cond_covar/{study}/want_genes.txt")])

    parsl.wait_for_current_tasks()

    # With specific tissue list, calculate conditioned covariances
    tiss_list= list(pd.read_table(tiss_file, header=None)[0])
    tiss_list_prefix = Path(tiss_file).stem

    tiss_cond_covar_res = []
    for cur_tiss in tiss_list:
        out_files = [File(f"../output/intermediate_data/condTWAS_cond_covar/{study}/{gene_name}__{cur_tiss}.txt.gz") for gene_name in cojo_out_gene_names]
        tiss_cond_covar_res.append(
                calc_tiss_cond_covar(cur_tiss, study, want_genes, 
                    inputs=pred_cond_snps_res,
                    outputs=out_files)
                )
    parsl.wait_for_current_tasks()

    # output of tiss_cond_covar_res is one gene, for each tissue in tiss_list, use these results to run spredixcan
    for tiss_covar_res in tiss_cond_covar_res:
        for gene_covar_res in tiss_covar_res.outputs:
            gene_covar_fname = gene_covar_res.result()

            gene_covar_search_string = Path(Path(gene_covar_fname).stem).stem

            search_res = re.search("(ENSG\d{1,}\.\d{1,})__(.*)", gene_covar_search_string)

            gene_name = search_res.group(1)
            tiss_name = search_res.group(2)
            output_file = f"../output/intermediate_data/condTWAS_spredixcan/{study}/spredixcan_igwas_gtexmashrv8_{gene_name}__PM__{tiss_name}.csv"
            spredixcan(gene_name, tiss_name, study, gene_covar_fname, output_file)

    # Some spredixcans can fail due to no tissue overlap. This is a workaround, since spredixcan doesn't return a DataFuture
    parsl.wait_for_current_tasks()

    # if splicing (no expression in name), run SMulTiXcanByFeature.py (ACAT),
    # otherwise, can run 1step acat in a python app.
    if "expression" in study:
        # def calc_1step_acat(results_dir, out_dir, tiss_list=""
        #         tissue_pattern="__PM__(.*)\.csv",
        #         gene_pattern="gtexmashrv8_(.*)__PM__", replace_pval=True,
        #         has_breast=True, outputs=[],):
        out_dir = f"../output/condTWAS_final/"
        final_file = os.path.join(out_dir, f"{study}_{tiss_list_prefix}.txt")
        combine_eqtl_and_calc_acat(
                study=study,
                results_dir=f"../output/intermediate_data/condTWAS_spredixcan/{study}/",
                out_dir=out_dir,
                tiss_list=tiss_file,
                has_breast=has_breast,
                outputs=[File(final_file)])
    else:
        sqtl_acat_single_gene_res = []
        for gene_parse_res, gene_name in zip(gene_parsed_res, cojo_out_gene_names):
            gene_parse_file = gene_parse_res.outputs[0]
            sqtl_acat_out_file = File(f"../output/intermediate_data/condTWAS_single_gene_mtiss_acat/{study}/{tiss_list_prefix}___ixcan_igwas_gtex_mashr_v8_ccn30_{gene_name}.txt")

            sqtl_acat_single_gene_res.append(
                    calc_mtiss_acat(gene_name=gene_name, study=study, tiss_list=tiss_file,
                                    inputs=[gene_parse_file],
                                    outputs=[sqtl_acat_out_file])

            )
        combine_sqtl_acat(study=study,
                results_dir=f"../output/intermediate_data/condTWAS_single_gene_mtiss_acat/{study}/",
                out_dir=f"../output/condTWAS_final/",
                tiss_list=tiss_file,
                has_breast=has_breast,
                inputs=sqtl_acat_single_gene_res,
                outputs=[File("f{study}_{tiss_list_prefix}.txt")])
    parsl.wait_for_current_tasks()
    return()


BRCA_INDEX_SNP_WINDOW = 2000000 # 2 MB
BRCA_CLUMP_R2 = 0.65

# Overall brca. BCAC and BCAC/UKB metal
# run_cojo("BCAC_Overall_BreastCancer_EUR", index_snp_window=BRCA_INDEX_SNP_WINDOW, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("BCAC_Overall_BreastCancer_EUR_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("meta_analysis_BCAC_UKB_ovr_brca", index_snp_window=BRCA_INDEX_SNP_WINDOW, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("meta_analysis_BCAC_UKB_ovr_brca_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
# 
# 
# # Erneg and erpos
# run_cojo("meta_analysis_BCAC_CIMBA_erneg_brca", index_snp_window=BRCA_INDEX_SNP_WINDOW, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("meta_analysis_BCAC_CIMBA_erneg_brca_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("BCAC_ERPOS_BreastCancer_EUR", index_snp_window=BRCA_INDEX_SNP_WINDOW, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("BCAC_ERPOS_BreastCancer_EUR_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
# 
# 
# # Erneg and erpos matching
# run_cojo("meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch", index_snp_window=BRCA_INDEX_SNP_WINDOW, clump_r2_imputed=BRCA_CLUMP_R2,
#         ermatch="ern")
# run_cojo("meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2, ermatch="ern")
# run_cojo("BCAC_ERPOS_BreastCancer_EUR_index_ermatch",
#         index_snp_window=BRCA_INDEX_SNP_WINDOW, clump_r2_imputed=BRCA_CLUMP_R2,
#         ermatch="erp")
# run_cojo("BCAC_ERPOS_BreastCancer_EUR_index_ermatch_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2, ermatch="erp")
# # 
# # # Erneg and erpos swapped
# run_cojo("meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap", index_snp_window=BRCA_INDEX_SNP_WINDOW, clump_r2_imputed=BRCA_CLUMP_R2,
#         ermatch="erp")
# # 
# run_cojo("meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2, ermatch="erp")
# # 
# run_cojo("BCAC_ERPOS_BreastCancer_EUR_index_erswap",
#         index_snp_window=BRCA_INDEX_SNP_WINDOW, clump_r2_imputed=BRCA_CLUMP_R2,
#         ermatch="ern")
# # 
# run_cojo("BCAC_ERPOS_BreastCancer_EUR_index_erswap_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2, ermatch="ern")

# Debug
# run_ctwas("meta_analysis_BCAC_UKB_ovr_brca")
# run_ctwas("BCAC_ERPOS_BreastCancer_EUR_index_erswap")

# BRCA intrinsic subtype
# run_cojo("intrinsic_subtype_1", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=False, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("intrinsic_subtype_2", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=False, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("intrinsic_subtype_3", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=False, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("intrinsic_subtype_4", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=False, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("intrinsic_subtype_5", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=False, clump_r2_imputed=BRCA_CLUMP_R2)

# run_cojo("intrinsic_subtype_1_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
# run_cojo("intrinsic_subtype_2_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
run_cojo("intrinsic_subtype_3_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
run_cojo("intrinsic_subtype_4_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
run_cojo("intrinsic_subtype_5_expression", index_snp_window=BRCA_INDEX_SNP_WINDOW, expression=True, clump_r2_imputed=BRCA_CLUMP_R2)
parsl.wait_for_current_tasks()
