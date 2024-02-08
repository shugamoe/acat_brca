#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import parsl
import subprocess
from subprocess import call
from parsl.app.app import bash_app, python_app, join_app
from parsl.dataflow.futures import AppFuture

from config import id_for_memo_File
from helpers import SetupCOJO
setup = SetupCOJO()

COND_TWAS_EXECUTORS=["gp"]

@bash_app(cache=True, executors=COND_TWAS_EXECUTORS)
def parse_gene_gwas(gene_file, gene_name, study, inputs=[], outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    import numpy as np
    import pandas as pd
    if not os.path.exists(os.path.dirname(outputs[0].filepath)):
        os.makedirs(os.path.dirname(outputs[0].filepath), exist_ok=True)
    if os.path.exists(outputs[0].filepath):
        return("echo 'Output exists. Remove it or delete it.'")

    output = outputs[0].filepath
    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_parsing.py \
    -gwas_file {gene_file} \
    -output_column_map SNP variant_id \
    -output_column_map refA effect_allele \
    -output_column_map non_effect_allele non_effect_allele \
    -output_column_map bC effect_size \
    -output_column_map bC_se standard_error \
    -output_column_map Chr chromosome --chromosome_format -output_column_map bp position \
    -output_column_map pC pvalue \
    -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
    -snp_reference_metadata ../../00_prep_inputs_for_metaxcan/input/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA  \
    -output {output}
    """
    return(bash_command)

@bash_app(cache=True, executors=COND_TWAS_EXECUTORS)
def make_pred_and_cond_snps_gene_gwas(gene_name, study, inputs=[], outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    import numpy as np
    import pandas as pd
    if not os.path.exists(os.path.dirname(outputs[0].filepath)):
        os.makedirs(os.path.dirname(outputs[0].filepath), exist_ok=True)
    if (os.path.exists(outputs[0].filepath) and os.path.exists(outputs[1].filepath)):
        return("echo 'Outputs exists. Remove it or delete it.'")

    gene_parse_file = inputs[0]
    out_pred_file = outputs[0].filepath
    out_cond_file = outputs[1].filepath
    bash_command = \
    f"""
    PRE_COJO_SNPS=../output/intermediate_data/study_snp_lists/{study}/{gene_name}.snplist.all.for.prediction
    POST_COJO_SNPS={gene_parse_file}
    OUT_PRED_FILE={out_pred_file}

    SEARCH_STRING=`echo $(cat ${{PRE_COJO_SNPS}} | tr '\n' '|' | sed 's/|/_|/g')$(head -1 ${{PRE_COJO_SNPS}})_`

    if test -f ${{OUT_PRED_FILE}}; then
        echo "${{OUT_PRED_FILE}} already written"
    else
        zcat ${{POST_COJO_SNPS}} | grep -P "${{SEARCH_STRING}}" | cut -f 2 > ${{OUT_PRED_FILE}}
        echo "Wrote ${{OUT_PRED_FILE}}"
    fi


    OUT_COND_FILE={out_cond_file}
    CLUMP_COND_COJO="../output/intermediate_data/clump_cond_snp_lists/{study}/{gene_name}.clumped.trimmed"
    COJO_MA=../output/intermediate_data/cojo_ma_files/{study}_gwas.ma
    SEARCH_STRING=`echo $(cat ${{CLUMP_COND_COJO}} | tr '\n' '|' | sed 's/|/ |/g')$(head -1 ${{CLUMP_COND_COJO}}) `
    
    if test -f ${{OUT_COND_FILE}}; then
        echo "${{OUT_COND_FILE}} already written"
    else
        cat ${{COJO_MA}} | grep -P "${{SEARCH_STRING}}" | awk '{{print $1"_"$3"_"$2"_b38"}}' > ${{OUT_COND_FILE}}
        echo "Wrote ${{OUT_COND_FILE}}"
    fi
    """
    return(bash_command)

@bash_app(cache=True, executors=COND_TWAS_EXECUTORS)
def calc_tiss_cond_covar(tiss_name, study, want_genes, inputs=[], outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    import numpy as np

    outdir = os.path.dirname(outputs[0].filepath)

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    if os.path.exists(outputs[0].filepath):
        return("echo 'Output exists. Remove it or delete it.'")
    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/conditional_covariance_for_model_group.py \
    -parquet_genotype_folder ../../00_prep_inputs_for_metaxcan/input/gtex_v8_parquet_eur_maf0.01_biallelic/ \
    -parquet_genotype_pattern "gtex_v8_eur_itm.chr(\d+).variants.parquet" \
    -group ../input/combine_phenotype_groups_multi_tiss.txt.gz \
    -model_db_group_key ../input/combine_eqtl_for_covar_calc.db \
    -model_db_group_value ../input/combine_sqtl_for_covar_calc.db \
    -parsimony 8 \
    -want_genes {want_genes} \
    -condition_info_dir ../output/intermediate_data/condTWAS_cond_covar/{study} \
    -pred_pattern {{want_gene}}.snps.pred.txt \
    -cond_pattern {{want_gene}}.snps.condition.txt \
    --individuals ../../02_acat_eqtl_sqtl/input/gtex_v8_by_tiss_individuals/{tiss_name}_individuals.txt \
    -tissue {tiss_name} \
    -output_dir {outdir}
    """

    bash_command_expression = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/groups_and_conditioned_covariance_for_model.py \
    --parquet_genotype_folder ../../00_prep_inputs_for_metaxcan/input/gtex_v8_parquet_eur_maf0.01_biallelic/ \
    --parquet_genotype_pattern "gtex_v8_eur_itm.chr(\d+).variants.parquet" \
    --model_db ../input/combine_eqtl_for_covar_calc.db \
    --parsimony 8 \
    --want_genes {want_genes} \
    --condition_info_dir ../output/intermediate_data/condTWAS_cond_covar/{study} \
    --pred_pattern {{want_gene}}.snps.pred.txt \
    --cond_pattern {{want_gene}}.snps.condition.txt \
    --individuals ../../02_acat_eqtl_sqtl/input/gtex_v8_by_tiss_individuals/{tiss_name}_individuals.txt \
    --tissue {tiss_name} \
    --output_dir {outdir} \
    --output {outdir}/{tiss_name}.sigma11.condition.numbers.tsv \
    --covar_mode
    """

    if "expression" in study:
        return(bash_command_expression)
    else:
        return(bash_command)

@python_app(cache=True, executors=COND_TWAS_EXECUTORS)
def make_want_genes(study, inputs=[], outputs=[],
                    stdout=parsl.AUTO_LOGNAME, 
                    stderr=parsl.AUTO_LOGNAME): 
        import os
        import re
        import pandas as pd

        cond_covar_dir = f"../output/intermediate_data/condTWAS_cond_covar/{study}/"
        if not os.path.exists(cond_covar_dir):
            os.makedirs(cond_covar_dir, exist_ok=True) 

        fpaths = [file_obj.filepath for file_obj in inputs]
        want_genes_l = [re.search("(ENSG\d{1,}\.\d{1,})", fpath).groups(1)[0] for fpath in fpaths]
        want_genes_df = pd.DataFrame({"gene": want_genes_l})
        want_genes = outputs[0].filepath
        want_genes_df.to_csv(want_genes, header=False, index=False)
        return(want_genes)



@bash_app(cache=True, executors=COND_TWAS_EXECUTORS)
def spredixcan(gene, tissue, study, covar_file, output_file, inputs=[], outputs=[],
               stdout=parsl.AUTO_LOGNAME, 
               stderr=parsl.AUTO_LOGNAME): 
    import os
    import numpy as np

    outdir = os.path.dirname(output_file)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    if os.path.exists(output_file):
        return("echo 'Output exists. Remove it or delete it.'")
    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/SPrediXcan.py \
    --gwas_file ../output/intermediate_data/condTWAS_parse/{study}/{gene}.txt.gz \
    --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
    --model_db_path ../../01_spredixcan_eqtl_sqtl/input/gtex_v8_sqtl_dbs_mashr/mashr_{tissue}.db \
    --covariance {covar_file} \
    --keep_non_rsid --additional_output --model_db_snp_key varID \
    --throw \
    --output_file {output_file}
    """

    bash_command_expression = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/SPrediXcan.py \
    --gwas_file ../output/intermediate_data/condTWAS_parse/{study}/{gene}.txt.gz \
    --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
    --model_db_path ../../01_spredixcan_eqtl_sqtl/input/gtex_v8_eqtl_dbs_mashr/mashr_{tissue}.db \
    --covariance {covar_file} \
    --keep_non_rsid --additional_output --model_db_snp_key varID \
    --throw \
    --output_file {output_file}
    """

    if "expression" in study:
        return(bash_command_expression)
    else:
        return(bash_command)





@python_app(cache=True, executors=COND_TWAS_EXECUTORS)
def combine_eqtl_and_calc_acat(study, results_dir,
                    out_dir,
                    tiss_list,
                    tissue_pattern="__PM__(.*)\.csv",
                    gene_pattern="gtexmashrv8_(.*)__PM__",
                    replace_pval=True,
                    has_breast=True,
                    outputs=[],
                    stdout=parsl.AUTO_LOGNAME, 
                    stderr=parsl.AUTO_LOGNAME): 

    import re
    import numpy as np
    import pandas as pd
    import copy
    from glob import glob

    if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    def get_gene_cojo_status(study):
        whitelist = pd.read_csv(f"../output/intermediate_data/cojo_input/{study}_whitelist.tsv", sep="\t", usecols=["group", "status"]
                )
        blacklist = pd.read_csv(f"../output/intermediate_data/cojo_input/{study}_blacklist.tsv", sep="\t", usecols=["group", "status"]                )
        gene_cojo_status = pd.concat([whitelist, blacklist])
        gene_cojo_status.loc[:, "cojo_status"] = gene_cojo_status['status']

        return(gene_cojo_status[["group", "cojo_status"]])

    def single_acat(pval_dict, has_breast):
        group_size = len(pval_dict["pval"])
        pval_arr = np.array(pval_dict["pval"], dtype=np.float64)
        pval_arr = pval_arr[np.logical_not(np.isnan(pval_arr))]
        available_results = pval_arr.shape[0]

        s = np.sum(np.tan((0.5 - pval_arr) * np.pi))
        acat = 0.5 - np.arctan(s) / np.pi
        if has_breast is True:
            try:
                breast_index = pval_dict["tissue"].index("Breast_Mammary_Tissue")
                breast_pvalue = pval_dict["pval"][breast_index]

                # Get an ACAT p-value calculated without breast tissue
                found_breast = True
                no_breast_pvals = copy.deepcopy(pval_dict["pval"])
                no_breast_pvals.pop(breast_index)

                no_breast_pval_arr = np.array(no_breast_pvals, dtype=np.float64)
                no_breast_pval_arr = no_breast_pval_arr[np.logical_not(np.isnan(no_breast_pval_arr))]

                s = np.sum(np.tan((0.5 - no_breast_pval_arr) * np.pi))
                no_breast_acat = 0.5 - np.arctan(s) / np.pi
            except ValueError:
                found_breast = False
                breast_pvalue = None
                no_breast_acat = acat


        if has_breast is True:
            return acat, no_breast_acat, breast_pvalue, group_size, available_results
        else:
            return acat, group_size, available_results

    rfiles = glob(os.path.join(results_dir, "*.csv"))
    pvals_dict = {}
    if tiss_list is not None:
        whitelist = list(pd.read_csv(tiss_list, header=None)[0])
    for gene_tiss_file in rfiles:
        tissue_regexp = re.compile(tissue_pattern)
        gene_regexp = re.compile(gene_pattern)
        tissue = tissue_regexp.search(os.path.basename(gene_tiss_file)).groups()[0]
        if tiss_list is not None:
            if tissue not in whitelist:
                continue
            else:
                pass
                # print("Using tissue: " + tissue)

        gene = gene_regexp.search(os.path.basename(gene_tiss_file)).groups()[0]

        df = pd.read_csv(gene_tiss_file, usecols=["gene", "pvalue"])
        df = df[df["gene"] == gene]
        if len(df) == 0:
            continue
        
        pval = df["pvalue"].tolist()[0]
        if replace_pval is True:
            if pval == 1:
                pval = pval - 1e-5

        if gene not in pvals_dict:
            pvals_dict[gene] = {"tissue": [tissue],
                                "pval": [pval]}
        elif tissue not in pvals_dict[gene]["tissue"]:
            pvals_dict[gene]["tissue"].append(tissue)
            pvals_dict[gene]["pval"].append(pval)
        else:
            print("What's going on?")
            print(tissue)
            print(pval)

    if has_breast is True:
        out_cols = ["group", "no_breast_acat", "acat", "breast_pvalue", "group_size", "available_results"]
    else:
        out_cols = ["group", "acat", "group_size", "available_results"]
    results = []
    for gene in pvals_dict:
        if has_breast is True:
            acat, no_breast_acat, breast_pvalue, group_size, available_results = single_acat(pvals_dict[gene], has_breast=has_breast)
            results.append((gene, acat, no_breast_acat, breast_pvalue, group_size, available_results))
        else:
            acat, group_size, available_results = single_acat(pvals_dict[gene], has_breast=has_breast)
            results.append((gene, acat, group_size, available_results))
    results = pd.DataFrame(results, columns=out_cols)

    merged = results.merge(get_gene_cojo_status(study), left_on="group", right_on="group", how="right")

    if has_breast is True:
        merged.loc[merged['acat'].isna(), ['acat', 'no_breast_acat', 'breast_pvalue']] = merged['cojo_status']
        merged = merged[["group", "acat", "no_breast_acat", "breast_pvalue", "group_size", "available_results"]]
    else:
        merged.loc[merged['acat'].isna(), ['acat']] = merged['cojo_status']
        merged = merged[["group", "acat", "group_size", "available_results"]]

    tiss_basename = os.path.basename(tiss_list)
    merged.to_csv(os.path.join(out_dir, f"{study}_{tiss_basename}"), index=False, sep="\t")
    return

@bash_app(cache=True, executors=COND_TWAS_EXECUTORS)
def calc_mtiss_acat(gene_name, study, tiss_list, inputs=[], outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    import numpy as np
    import pandas as pd
    if not os.path.exists(os.path.dirname(outputs[0].filepath)):
        os.makedirs(os.path.dirname(outputs[0].filepath), exist_ok=True)
    if (os.path.exists(outputs[0].filepath)):
        return("echo 'Output exists. Remove it or delete it.'")

    gene_parse_file = inputs[0]
    out_file = outputs[0]
    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/SMulTiXcanByFeature.py \
    --gwas_file {gene_parse_file} \
    --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore  \
    --grouping ../../02_acat_eqtl_sqtl/input/combine_phenotype_groups.txt.gz GTEx_sQTL \
    --model_db_path ../../02_acat_eqtl_sqtl/input/combine_sqtl.db \
    --keep_non_rsid --model_db_snp_key varID \
    --acat \
    --tiss_list {tiss_list} \
    --covariance ../../02_acat_eqtl_sqtl/input/combine_cross_intron_covar.txt.gz \
    --associations_folder ../output/intermediate_data/condTWAS_spredixcan/{study}/ \
    --associations_file_pattern spredixcan_igwas_gtexmashrv8_{gene_name}* \
    --associations_tissue_pattern "__PM__(.*)\.csv" \
    --single_gene "{gene_name}" \
    --cutoff_condition_number 30 \
    --verbosity 10 \
    --output {out_file}
    """
    return(bash_command)




@python_app(cache=True, executors=COND_TWAS_EXECUTORS)
# def combine_eqtl_and_calc_acat(study, results_dir,
#                     out_dir,
#                     tiss_list,
def combine_sqtl_acat(study, results_dir, out_dir, tiss_list,
        has_breast,
        inputs=[],
        outputs=[],
        stdout=parsl.AUTO_LOGNAME, 
        stderr=parsl.AUTO_LOGNAME): 
    import os
    import pandas as pd
    from glob import glob
    from pathlib import Path
    
    def get_gene_cojo_status(study):
        whitelist = pd.read_csv(f"../output/intermediate_data/cojo_input/{study}_whitelist.tsv", sep="\t", usecols=["group", "status"])
        blacklist = pd.read_csv(f"../output/intermediate_data/cojo_input/{study}_blacklist.tsv", sep="\t", usecols=["group", "status"])
        gene_cojo_status = pd.concat([whitelist, blacklist])
        gene_cojo_status.loc[:, "cojo_status"] = gene_cojo_status['status']

        return(gene_cojo_status[["group", "cojo_status"]])

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    tiss_list_name = Path(tiss_list).stem

    result_paths = glob(os.path.join(results_dir, f"{tiss_list_name}___*.txt"))
    rfiles = glob(os.path.join(results_dir, "*.csv"))
    print(f"Using tissue list: {tiss_list_name}")
    print(f"{len(result_paths)} results found")

    made_base = False
    for rpath in result_paths:
        cur_res = pd.read_csv(rpath, sep="\t")
        if made_base is False:
            base = cur_res
            made_base = True
        else:
            base = base.append(cur_res)

    write_df = base.merge(get_gene_cojo_status(study), left_on="group", right_on="group", how="left")

# No longer joining to original results
# write_df = write_df.merge(pre_cojo_res, left_on="group", right_on="group", suffixes=('_conditioned', '_original'))

    if has_breast is True:
        write_df.loc[write_df['mtiss_acat'].isna(), ['mtiss_acat', 'breast_acat']] = write_df['cojo_status']
        write_df.loc[write_df['mtiss_acat'] == 'has_original_index_snps', ['mtiss_acat', 'breast_acat']] = 'ref_sample_conditional_freq_mismatch'
    else:
        write_df.loc[write_df['mtiss_acat'].isna(), ['mtiss_acat']] = write_df['cojo_status']
        write_df.loc[write_df['mtiss_acat'] == 'has_original_index_snps', ['mtiss_acat']] = 'ref_sample_conditional_freq_mismatch'

    write_df.to_csv(os.path.join(out_dir, f"{study}_{tiss_list_name}.txt"), index=False, sep="\t")
    return()
