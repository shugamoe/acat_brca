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

@bash_app(cache=True)
def make_study_snp_lists(study, gene_whitelist_colname, inputs=[], outputs=[],
                        expression=False,
                        stdout=parsl.AUTO_LOGNAME, 
                        stderr=parsl.AUTO_LOGNAME): 
    import os
    if not os.path.exists(setup.STUDY_SNP_LIST_DIR_KEY[study]):
        os.makedirs(setup.STUDY_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if os.path.exists(outputs[0].filepath): 
        return("echo 'Output exists. Remove it or delete it.'")
    cojo_ma_file = inputs[0]
    gene_whitelist  = inputs[1]
    output = outputs[0].filepath

        # python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/groups_and_conditioned_covariance_for_model.py \
    if expression:
        bash_command = \
        f"""
        python3 /gpfs/data/huo-lab/jmcclellan/software/summary-gwas-imputation/src/groups_and_conditioned_covariance_for_model.py \
        --gwas_file {setup.GWAS_KEY[study]} \
        --model_db {setup.MODEL_EQTL_DB} \
        --output_dir {setup.STUDY_SNP_LIST_DIR_KEY[study]} \
        --output {output} \
        --parsimony 7 \
        --get_og_for_imputed {cojo_ma_file} \
        --gene_whitelist {gene_whitelist} "{gene_whitelist_colname}"
        """
    else:
        # python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/SMulTi_gene_snps.py \
        bash_command = \
        f"""
        echo $PATH
        which python3
        python3 /gpfs/data/huo-lab/jmcclellan/software/MetaXcan/software/SMulTi_gene_snps.py \
        --gwas_file {setup.GWAS_KEY[study]} \
        --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
        --grouping {setup.GROUPING} GTEx_sQTL \
        --model_db_path {setup.MODEL_SQTL_DB} \
        --keep_non_rsid --model_db_snp_key varID \
        --associations {setup.MTISS_TWAS_KEY[study]} SPrediXcan \
        --verbosity 10 \
        --get_og_for_imputed {cojo_ma_file} \
        --snp_list_prefix {setup.STUDY_SNP_LIST_DIR_KEY[study]} \
        --output {output} \
        --gene_whitelist {gene_whitelist} "{gene_whitelist_colname}"
        """
    return(bash_command) 

@python_app(cache=True)
def make_cojo_ma_file(study, outputs=[],
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
    pre_output = os.path.join(os.path.dirname(output), "not_final_" + os.path.basename(output))
    bash_command = \
    f"""
    zcat {setup.GWAS_KEY[study]} | awk 'BEGIN {{FS="\t"}}; {{print $3"_"$4" "$5" "$6" "$8" "$12" "$13" "$11" "$10" {setup.SAMPLE_N_KEY[study]}"}}' > {pre_output}
    """
    proc = subprocess.run(bash_command, shell=True)
    if proc.returncode == 0:
        try:
            pre_cojo_ma_names = ["chromosome_position", "effect_allele", "non_effect_allele", "frequency", "effect_size", "standard_error", "pvalue", "zscore", "N"]
            pre_cojo_ma_df = pd.read_csv(pre_output, header=0, names=pre_cojo_ma_names, delim_whitespace=True)

            pre_cojo_ma_df['effect_size_alt'] = pre_cojo_ma_df['zscore'] / np.sqrt(setup.SAMPLE_N_KEY[study])
            pre_cojo_ma_df['N'][pre_cojo_ma_df['N'].isnull()] = setup.SAMPLE_N_KEY[study]

            pre_cojo_ma_df['effect_size'][pre_cojo_ma_df['effect_size'].isnull()] = pre_cojo_ma_df['effect_size_alt'][pre_cojo_ma_df['effect_size'].isnull()]
            pre_cojo_ma_df['standard_error'][pre_cojo_ma_df['standard_error'].isnull()] = 1 / np.sqrt(setup.SAMPLE_N_KEY[study]) 

            cojo_ma_names = ["chromosome_position", "effect_allele", "non_effect_allele", "frequency", "effect_size", "standard_error", "pvalue", "N"]
            pre_cojo_ma_df[cojo_ma_names].to_csv(output, sep=" ", header=True, index=False)
            os.remove(pre_output)
        except FileNotFoundError:
            print("Shouldn't have no file.")
    else:
        raise Exception(f"shell command failed? {proc}")
    return

@python_app(cache=True)
def make_cond_snp_lists(study, inputs=[], outputs=[], bp_range=5000000, imp_check_range=500000, ermatch=None):
                       # stdout=parsl.AUTO_LOGNAME, 
                       # stderr=parsl.AUTO_LOGNAME): 
    def get_gencode(path=setup.GENCODE):
        import pandas as pd
        gencode = pd.read_csv(path, sep="\t", usecols=["chromosome", "start_location", "end_location", "gene_id", "gene_name"])
        gencode["chromosome"] = gencode["chromosome"].str.extract(r"chr(\d{1,2})")
        gencode = gencode.dropna() # Leave out X, Y, M chromosomes (NaN)
        gencode["chr_num"] = gencode["chromosome"].astype(int) # int64 conversion
        return gencode[["chr_num", "start_location", "end_location", "gene_id", "gene_name"]]
    
    def get_snplist(path=setup.REPORTED_SNP_LIST, ermatch=None):
        import numpy as np
        if path == setup.REPORTED_SNP_LIST:
            snplist = pd.read_csv(path, sep="\t")
            if ermatch == "erp":
                snplist = snplist[snplist['ER+ p'] < 5e-8]
            elif ermatch == "ern":
                snplist = snplist[snplist['ER- p'] < 5e-8]
            else:
                pass
            snplist['chr_num'] = snplist['Chromosome']
            snplist['hg38_SNP'] = snplist.agg('chr{0[Chromosome]}_{0[hg38 position]}'.format, axis=1)
            snplist['Position (hg38)'] = snplist['hg38 position']

        # hg38_SNP is of the form chr#_<pos>
        return snplist[["chr_num", "hg38_SNP", "Position (hg38)"]]

    import os
    import pandas as pd
    from scipy import stats
    import numpy as np

    if not os.path.exists(setup.COND_SNP_LIST_DIR_KEY[study]):
        os.makedirs(setup.COND_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if not os.path.exists(setup.PROXY_COND_SNP_LIST_DIR_KEY[study]):
        os.makedirs(setup.PROXY_COND_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if os.path.exists(outputs[0].filepath):
        print("Output exists. Remove it or delete it.")
        return

    gencode = get_gencode()
    snplist = get_snplist(ermatch=ermatch)

    study_genes = pd.read_csv(inputs[0].filepath, sep="\t")
    cojo_ma_df = pd.read_csv(inputs[1].filepath, sep=" ", usecols=['chromosome_position', 'effect_size', 'standard_error'])
    cojo_ma_df['position'] = cojo_ma_df['chromosome_position'].str.split("_", expand=True)[1].astype(int)
    cojo_ma_df['chr_num'] = cojo_ma_df["chromosome_position"].str.extract(r"(\d{1,2})").astype(int)
    
    cond_snps_tracker = []
    for gene in study_genes['group']:
        # snplist.all.for.prediction if gene has imputed and original model snps
        # only .snplist if gene only has original model SNPs
        primary_path = "{}".format(os.path.join(setup.STUDY_SNP_LIST_DIR_KEY[study], f"{gene}.snplist.all.for.prediction"
            ))
        secondary_path = "{}".format(os.path.join( setup.STUDY_SNP_LIST_DIR_KEY[study], f"{gene}.snplist"
            ))
        if os.path.exists(primary_path):
            use_path = primary_path
        else:
            use_path = secondary_path

        gene_model_snps = pd.read_csv(use_path, header=None, names=["model_snp"])
        nearby_og_snps_path = 'NA'
        cond_snps_path = 'NA'
        gene_info = gencode[gencode["gene_id"] == gene] # gene_id is ENSG*, gene_name is "WASH7P" sort of stuff
        if len(gene_info) == 0:
            print(f"Can't find this gene in the gencode: {setup.GENCODE}")
            status = "not_found"
            cond_snps_tracker.append((gene, cond_snps_path, None, status, nearby_og_snps_path))
            continue
        elif len(gene_info) > 1:
            print(f"Multiple gencode matches for {gene}, probably repeats, just using the first row.")
            print(gene_info)
            gene_info = gene_info.loc[gene_info.index == 0]
            # status = "multiple_matches" # No index SNPs to condition on
            # cond_snps_tracker.append((gene, cond_snps_path, None, status, nearby_og_snps_path))
            continue

        gene_chr_num = int(gene_info['chr_num'])
        start_loc = gene_info['start_location'] - bp_range
        end_loc = gene_info['end_location'] + bp_range
        cond_snps = snplist[snplist["Position (hg38)"].between(int(start_loc), int(end_loc))]
        cond_snps = cond_snps[cond_snps['chr_num'] == gene_chr_num]
        cond_snps = cond_snps.drop_duplicates(subset=["hg38_SNP"]) # Remove dups if there for some reason

        if len(cond_snps) == 0:
            print(f"{gene} does't have any reported SNPs within {bp_range} base pairs")
            status = "no_index" # No index SNPs to condition on
            cond_snps_tracker.append((gene, cond_snps_path, gene_chr_num, status, nearby_og_snps_path))
            continue

        cond_snps = cond_snps.merge(cojo_ma_df[['chromosome_position',
            'position', 'effect_size', 'standard_error']], left_on='hg38_SNP',
            right_on='chromosome_position', how='inner')
        cond_snps_imp = cond_snps[cond_snps['effect_size'].isna()] # Imputed SNPS
        cond_snps_og = cond_snps[cond_snps['effect_size'].notna()] # Typed, original SNPs

        # If any of our model snps for the gene overlap with cond_snps_og
        # (original snps used as initial list of conditional snps), save these
        # overlapping SNPs. We will want to 0 them out later after we get COJO
        # results
        gene_model_cond_snps_olap = gene_model_snps.merge(cond_snps_og, left_on='model_snp', 
                right_on='hg38_SNP', how='inner')
        if gene_model_cond_snps_olap.shape[0] > 0:
            model_cond_olap_path = f"{setup.COND_SNP_LIST_DIR_KEY[study]}{gene}.snplist.model.cond.overlap"
            gene_model_cond_snps_olap[["model_snp"]].to_csv(model_cond_olap_path, sep=" ", header=True, index=False)

        if cond_snps_imp.shape[0] > 0:
            # For inspection later, save the original typed SNPs if we're going to add in proxies later for imputed
            cond_snps_path = f"{setup.COND_SNP_LIST_DIR_KEY[study]}{gene}.og.no.proxy.snplist"
            cond_snps_og['zscore'] = cond_snps_og.effect_size / cond_snps_og.standard_error
            cond_snps_og['P'] = 2 * stats.norm.sf(np.abs(cond_snps_og.zscore))
            cond_snps_og = cond_snps_og.sort_values(by='P', ascending=False) # Lowest Pvalue on bottom
            cond_snps_og[["hg38_SNP", "P"]].to_csv(cond_snps_path, sep=" ", header=True, index=False)

            # Make a list of nearby OG snps
            all_nearby_og_snps = None
            cojo_ma_df_slice = cojo_ma_df[cojo_ma_df['chr_num'] == gene_chr_num]
            cojo_ma_df_slice = cojo_ma_df_slice.dropna() # Want original typed genes only
            for imp_snp in cond_snps_imp.itertuples():
                imp_snp_nearby_og_snps = cojo_ma_df_slice[cojo_ma_df_slice['position'].between(imp_snp.position - imp_check_range - 1, imp_snp.position + imp_check_range + 1)]
                if all_nearby_og_snps is None:
                    all_nearby_og_snps = imp_snp_nearby_og_snps
                else:
                    all_nearby_og_snps = all_nearby_og_snps.append(imp_snp_nearby_og_snps)
                    all_nearby_og_snps = all_nearby_og_snps.drop_duplicates()

            # Don't want to use a model snp for the gene as a potential proxy
            # (low probability, but need to be sure)
            all_nearby_og_snps = all_nearby_og_snps[np.logical_not(
                all_nearby_og_snps['chromosome_position'].isin(gene_model_cond_snps_olap['model_snp']))]

            # Write the list of nearby OG snps
            nearby_og_snps_path = f"{setup.PROXY_COND_SNP_LIST_DIR_KEY[study]}{gene}.nearbyog.snplist"
            all_nearby_og_snps['zscore'] = all_nearby_og_snps.effect_size / all_nearby_og_snps.standard_error
            all_nearby_og_snps['P'] = 2 * stats.norm.sf(np.abs(all_nearby_og_snps.zscore))
            all_nearby_og_snps = all_nearby_og_snps.sort_values(by='P', ascending=False) # Lowest Pvalue on bottom
            all_nearby_og_snps[["chromosome_position", "P"]].to_csv(nearby_og_snps_path, sep=" ", header=True, index=False)

            # Write the list of imputed SNPs
            cond_snps_imp_path = f"{setup.PROXY_COND_SNP_LIST_DIR_KEY[study]}{gene}.imputed.snplist"
            cond_snps_imp['P'] = 0.00000001
            cond_snps_imp[["chromosome_position", "P"]].to_csv(cond_snps_imp_path, sep=" ", header=True, index=False)

            bfile = setup.COJO_BFILE_PATTERN.format(chr_num=gene_chr_num)
            bash_command = f"""
            plink \
              --bfile {bfile} \
              --clump {cond_snps_imp_path},{nearby_og_snps_path} \
              --clump-r2 .7 \
              --clump-kb 500 \
              --clump-p1 0.0000001 \
              --clump-p2 0.000001 \
              --memory 16000 \
              --threads 1 \
              --clump-snp-field chromosome_position \
              --clump-field P \
              --clump-index-first \
              --clump-best \
              --out {setup.PROXY_COND_SNP_LIST_DIR_KEY[study]}{gene}
            """
            proc = subprocess.run(bash_command, shell=True)
            if proc.returncode == 0:
                try:
                    og_proxy_snps = pd.read_csv(f"{setup.PROXY_COND_SNP_LIST_DIR_KEY[study]}{gene}.clumped.best", delim_whitespace=True)
                    og_proxy_snps = og_proxy_snps[['PSNP']].astype(str)
                    more_cond_snps_og = og_proxy_snps.merge(cojo_ma_df[['chromosome_position',
                        'position', 'effect_size', 'standard_error']], left_on='PSNP',
                        right_on='chromosome_position', how='inner')
                    more_cond_snps_og = more_cond_snps_og[['chromosome_position', 'position', 'effect_size', 'standard_error']]
                    num_og = len(cond_snps_og)
                    cond_snps_og = pd.concat([cond_snps_og, more_cond_snps_og])
                    num_og_after = len(cond_snps_og)
                    cond_snps_og = cond_snps_og.drop_duplicates()
                    has_imp = True
                except FileNotFoundError: # Sometimes find no clump results
                    pass
                    has_imp = False
            else:
                raise Exception(f"shell command failed? {proc}")
        else:
            has_imp = False


        if len(cond_snps_og) == 0:
            print(f"{gene} does't have any original reported SNPs within {bp_range} base pairs")
            gene_start = int(gene_info['start_location'])
            gene_end = int(gene_info['end_location'])
            min_dist_from_start = abs(cond_snps["Position (hg38)"] - gene_start).min()
            min_dist_from_end = abs(cond_snps["Position (hg38)"] - gene_end).min()
            
            status = "index_snps_imputed" # All index SNPs are imputed in original data, can't run COJO
            if (min_dist_from_start > imp_check_range) or (min_dist_from_end > imp_check_range):
                status += "_closest_over_500_kb_from_gene" 
            cond_snps['P'] = 0.00001

            # # Make a list of nearby OG snps
            # all_nearby_og_snps = None
            # cojo_ma_df_slice = cojo_ma_df[cojo_ma_df.chr_num == gene_chr_num]
            # for imp_snp in cond_snps.itertuples():
            #     index_in_ma = cojo_ma_df_slice.index[cojo_ma_df_slice.chromosome_position == imp_snp.hg38_SNP]
            #     imp_snp_nearby_og_snps = cojo_ma_df_slice.loc[int(index_in_ma.values) - og_check_window:int(index_in_ma.values) + og_check_window]
            #     if all_nearby_og_snps is None:
            #         all_nearby_og_snps = imp_snp_nearby_og_snps
            #     else:
            #         all_nearby_og_snps = all_nearby_og_snps.append(imp_snp_nearby_og_snps)
            #         all_nearby_og_snps = all_nearby_og_snps.drop_duplicates()

            # nearby_og_snps_path = f"{setup.COND_SNP_LIST_DIR_KEY[study]}{gene}.nearbyog.snplist"
            # cond_snps_path = f"{setup.COND_SNP_LIST_DIR_KEY[study]}{gene}.imputed.snplist"
            # cond_snps[["hg38_SNP", "P"]].to_csv(cond_snps_path, sep=" ", header=True, index=False)
            # all_nearby_og_snps['zscore'] = all_nearby_og_snps.effect_size / all_nearby_og_snps.standard_error
            # all_nearby_og_snps['P'] = 2 * stats.norm.sf(np.abs(all_nearby_og_snps.zscore))
            # all_nearby_og_snps = all_nearby_og_snps.sort_values(by='P', ascending=False) # Lowest Pvalue on bottom
            # all_nearby_og_snps[["hg38_SNP", "P"]].to_csv(nearby_og_snps_path, sep=" ", header=True, index=False)
        else:
            status = "has_original_index_snps"
            if has_imp is True:
                status += "_with_proxies"
            cond_snps_og['zscore'] = cond_snps_og.effect_size / cond_snps_og.standard_error
            cond_snps_og['P'] = 2 * stats.norm.sf(np.abs(cond_snps_og.zscore))
            cond_snps_og = cond_snps_og.sort_values(by='P', ascending=False) # Lowest Pvalue on bottom
            cond_snps_path = f"{setup.COND_SNP_LIST_DIR_KEY[study]}{gene}.snplist"
            cond_snps_og['hg38_SNP'] = cond_snps_og['chromosome_position']
            cond_snps_og[["hg38_SNP", "P"]].to_csv(cond_snps_path, sep=" ", header=True, index=False)
        cond_snps_tracker.append((gene, cond_snps_path, gene_chr_num, status, nearby_og_snps_path))

    tracker_df = pd.DataFrame(cond_snps_tracker, columns=["group", "cond_snp_list_path", "chr_num", "status", "nearby_og_snps_path"])
    tracker_df.to_csv(outputs[0].filepath, sep="\t", index=False)
    return

@python_app(cache=True)
def make_cojo_run_data(study, study_list, cond_genes_list, gene_col="group", sep="\t", inputs=[], outputs=[]):
    import pandas as pd
    if not os.path.exists(setup.COJO_IN_DIR_KEY[study]):
        os.makedirs(setup.COJO_IN_DIR_KEY[study], exist_ok=True)

    clump_cond_snp_lists = inputs
    whitelist = pd.read_csv(study_list, sep=sep) # group
    cond_genes_df = pd.read_csv(cond_genes_list, sep=sep)
    cond_genes_df['clump_cond_snp_list_path'] = clump_cond_snp_lists
    whitelist = whitelist.merge(cond_genes_df, how='left') # cond_snp_list_path chr_num

    blacklist = whitelist[whitelist['clump_cond_snp_list_path'].isna()] # If no conditional SNPs then can't run COJO
    whitelist = whitelist[whitelist['clump_cond_snp_list_path'].notna()] 

    blist_path, wlist_path = outputs 
    blacklist.to_csv(blist_path, index=False, sep="\t")
    whitelist.to_csv(wlist_path, index=False, sep="\t")
    return

@bash_app(cache=True)
def make_clump_cond_snp_list(study, chr_num, plink_output, cond_snp_list,
                      clump_snp_field="hg38_SNP", 
                      clump_field="P",
                      clump_kb=10000,
                      clump_r2=.85,
                      clump_p1=1, # Index variants are chosen greedily, this helps avoid "no clumps made" results
                      cojo_ma_path=None,
                      collin_prevent_mode="plink", # How collinearity issues in COJO --cojo-cond is avoided.
                      outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    if not os.path.exists(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study]):
        os.makedirs(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if os.path.exists(outputs[0].filepath):
        return("echo 'Output exists. Remove it or delete it.'")
    bfile = setup.COJO_BFILE_PATTERN.format(chr_num=chr_num)

    final_output = outputs[0].filepath
    out_prefix, _ = os.path.splitext(plink_output)
    plink_bash_command = \
    f"""
    plink \
     --bfile {bfile} \
     --clump {cond_snp_list} \
     --clump-r2 {clump_r2} \
     --clump-kb {clump_kb} \
     --clump-p1 {clump_p1} \
     --memory 16000 \
     --threads 1 \
     --clump-snp-field {clump_snp_field} \
     --clump-field {clump_field} \
     --out {out_prefix}

    if test -f "{plink_output}"; then
      tail -n +2 {plink_output} | tr -s ' ' | cut -d ' ' -f 4 | grep 'chr' > {final_output}
    else
      echo "No plink output? {plink_output}" # Use the most significant old conditioned SNP instead I guess
      cat {cond_snp_list} | tail -n +2 | cut -d ' ' -f 1 | tail -1  > {final_output} 
    fi
    """

    cojo_ma_path = cojo_ma_path.filepath
    gcta_cond_snp_list= cond_snp_list.filepath.replace("snplist", "gcta_snplist")
    gcta_out = out_prefix + ".jma.cojo"
    gcta_final_out = out_prefix + ".clumped.trimmed"

    gcta_bash_command = \
    f"""
    cat {cond_snp_list} | cut -d' ' -f 1 | grep -v 'hg38_SNP' > {gcta_cond_snp_list}


    gcta-1.94.1 \
     --bfile {bfile} \
     --cojo-file {cojo_ma_path} \
     --extract {gcta_cond_snp_list} \
     --cojo-slct \
     --out {out_prefix}

    if test -f "{gcta_out}"; then
      cat {gcta_out} | cut -f 2 | grep -v 'SNP' > {gcta_final_out}
    else
      echo "No SNPs selected. Choosing Top 1"
      cat {gcta_cond_snp_list} | head -1 > {gcta_final_out}
    fi
    """
    if collin_prevent_mode == "gcta":
        return(gcta_bash_command)
    elif collin_prevent_mode == "plink":
        return(plink_bash_command)

@bash_app(cache=True)
def make_clump_og_proxy_cond_snp_list(study, chr_num, plink_output, cond_snp_list, nearby_og_snp_list,
                      clump_snp_field="hg38_SNP", 
                      clump_field="P",
                      clump_kb=10000,
                      clump_r2=.85,
                      clump_p1=1, # Index variants are chosen greedily, this helps avoid "no clumps made" results
                      clump_p2=1, # We don't care about the proxy's p-value
                      outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    if not os.path.exists(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study]):
        os.makedirs(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if os.path.exists(outputs[0].filepath):
        return("echo 'Output exists. Remove it or delete it.'")
    bfile = setup.COJO_BFILE_PATTERN.format(chr_num=chr_num)

    final_output = outputs[0].filepath
    out_prefix, _ = os.path.splitext(plink_output)
    bash_command = \
    f"""
    plink \
     --bfile {bfile} \
     --clump {cond_snp_list},{nearby_og_snp_list} \
     --clump-r2 {clump_r2} \
     --clump-kb {clump_kb} \
     --clump-p1 {clump_p1} \
     --clump-p2 {clump_p2} \
     --memory 16000 \
     --threads 1 \
     --clump-snp-field {clump_snp_field} \
     --clump-field {clump_field} \
     --clump-replicate \
     --clump-index-first \
     --clump-best \
     --out {out_prefix}

    if test -f "{plink_output}"; then
      tail -n +2 {plink_output} | tr -s ' ' | cut -d ' ' -f 2 | grep 'chr' > {final_output}
    else
      echo "No plink output? {plink_output}" # Use the most significant old conditioned SNP instead I guess
      cat {cond_snp_list} | tail -n +2 | cut -d ' ' -f 1 | tail -1  > {final_output} 
    fi
    """
    return(bash_command)

@python_app(cache=False)
def add_non_effect_allele(study, cojo_ma_path, cojo_res_sep="\t", cojo_ma_sep=" ", inputs=[]):
    import numpy as np
    import pandas as pd 
    from glob import glob
    import os
    import re
    cojo_ma_df = pd.read_csv(cojo_ma_path, sep=cojo_ma_sep)

    cojo_oput_path = inputs[0].filepath
    oput_df = pd.read_csv(cojo_oput_path, sep=cojo_res_sep)
    gene = re.split(".cma.cojo", os.path.basename(cojo_oput_path))[0] # Could be a "protein" or "sequence id" maybe soon too

    # Put in non effect allele if not present
    if 'non_effect_allele' not in oput_df:
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
        oput_df.to_csv(new_out, sep=cojo_res_sep, index=False)
    return

@python_app(cache=True)
def make_alt_snp_lists(study, na_snp_window, whitelist_path, cojo_ma_path, cojo_ma_sep=" ", cojo_res_sep="\t"):
    import pandas as pd
    from glob import glob
    import os
    cojo_ma_df = pd.read_csv(cojo_ma_path, sep=cojo_ma_sep, usecols=['chromosome_position', 'effect_size'])
    cojo_ma_df = cojo_ma_df.dropna()
    cojo_ma_df.reset_index(drop=True, inplace=True)
    cojo_ma_df['chr_num'] = cojo_ma_df["chromosome_position"].str.extract(r"(\d{1,2})").astype(int)
    cojo_ma_df = cojo_ma_df[['chromosome_position', 'chr_num']]

    na_results = glob(f"{setup.COJO_OUT_DIR_KEY[study]}*.cma.allna")
    whitelist = pd.read_csv(whitelist_path, sep="\t")
    for nar in na_results:
        fname_only = os.path.basename(nar) 
        gene = fname_only.split(".cma.allna")[0]
        na_df = pd.read_csv(nar, sep=cojo_res_sep, usecols=["SNP", "Chr"])
        cojo_ma_df_slice = cojo_ma_df[cojo_ma_df.chr_num == na_df['Chr'][0]]
        new_snplist = None
        for na_snp in na_df.itertuples():
            index_in_ma = cojo_ma_df_slice.index[cojo_ma_df_slice.chromosome_position == na_snp.SNP]
            na_snp_nearby_og_snps = cojo_ma_df_slice.loc[int(index_in_ma.values) - na_snp_window:int(index_in_ma.values) + na_snp_window]
            if new_snplist is None:
                new_snplist = na_snp_nearby_og_snps
            else:
                new_snplist = new_snplist.append(na_snp_nearby_og_snps)
                new_snplist = new_snplist.drop_duplicates()
        ofile = f"{setup.STUDY_SNP_LIST_DIR_KEY[study]}{gene}.snplist.nafix"
        new_snplist[['chromosome_position']].to_csv(ofile, index=False, header=False)
        whitelist.loc[whitelist.group == gene, 'na_fix_snp_list_path'] = ofile
        whitelist.loc[whitelist.group == gene, 'status'] = 'has_original_index_snps_na_fix' 
    written = whitelist.to_csv(whitelist_path, index=False, sep="\t")
    return(written)





