# Makes the tables as specified in ../input/tablesjulian_format.docx
# Tons of pair-wise copy pasted (eqtl/sqtl) code.
library(purrr)

STUDIES <- c("meta_analysis_BCAC_CIMBA_erneg_brca", "meta_analysis_BCAC_UKB_ovr_brca", "BCAC_Overall_BreastCancer_EUR",
"BCAC_ERPOS_BreastCancer_EUR")
CTWAS_STUDIES <- c("meta_analysis_BCAC_UKB_ovr_brca", paste0("intrinsic_subtype_", 1:5))
# CTWAS_STUDIES <- c()

ER_MATCH_STUDIES <- c("meta_analysis_BCAC_CIMBA_erneg_brca", "BCAC_ERPOS_BreastCancer_EUR")


make_tables <- function(study_name,
                 table_names = c("supp_table1__{study_name}.tsv", "table1_twas_unrep__{study_name}.tsv",
                                 "table2_cojo_sig__{study_name}.tsv",
                                 "table3_cand_causal__{study_name}.tsv",
                                 "supp_table2or3__{study_name}.tsv"),
                 tiss_filter = "../../02_acat_eqtl_sqtl/input/tissue_lists/11tiss.txt",
                 gcode_convert = "../output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps_2023Jun02.tsv",
                 odir = "../output/drilldown_tables/",
                 return_raw=F,
                 filter_sig=T,
                 test=F
                 ){
  require(tidyverse)
  require(data.table)
  require(glue)
  require(purrr)

  tiss_list_basename_no_ext <- basename(tiss_filter) %>%
    str_replace_all("\\.txt", "")

  print(glue("{study_name}"))

  if (!dir.exists(odir)){
    dir.create(odir, recursive=T)
  }

  # chromosome	start_location.v26	end_location.v26	feature_type.v26	strand	gene_id.v26	gene_name.v26	gene_type.v26	gene_id_no_decimal	source	feature_type.v40	start_location.v40	end_location.v40	gene_id.v40	gene_type.v40	gene_name.v40	gencode_conversion_success	start_location	end_location	Gene Loci	Index SNP	Index SNP Loci	distance_to_snp	distance_to_start	distance_to_stop	Distance from Index SNP(kb)	Distance from Gene Start(kb)	Distance from Gene End(kb)	brca_twas_reported
  gcode_df <- fread(gcode_convert) %>%
    mutate(genomic_coord = glue("{chromosome}:{start_location.v40}-{end_location.v40}")) %>%
    mutate(genomic_coord = str_replace(genomic_coord, "chr", "")) %>%
    mutate(`genomic_coord.v26` = glue("{chromosome}:{start_location.v26}-{end_location.v26}")) %>%
    mutate(`genomic_coord.v26` = str_replace(`genomic_coord.v26`, "chr", ""))

  # Gather <study_name> results, both eqtl and sqtl
  # Can be significant either from 11tiss or breast tissue bonferonni significance
  # group   acat    breast_pvalue   group_size  available_results
  acat_eqtl <- fread(glue("../../02_acat_eqtl_sqtl/output/acat_eqtl/11tiss__{study_name}_eqtl_acat_results.tsv"))

  thresh_eqtl <- .05 / (acat_eqtl %>% filter(!is.na(acat)) %>% nrow())
  thresh_eqtl_breast <- .05 / (acat_eqtl %>% filter(!is.na(breast_pvalue)) %>% nrow())

  acat_eqtl <- acat_eqtl %>%
    mutate(sig_eqtl = ifelse(acat < thresh_eqtl, T, F),
           sig_eqtl_breast = ifelse(breast_pvalue < thresh_eqtl_breast, T, F),
           sig_eqtl_any = ifelse((sig_eqtl == TRUE | sig_eqtl_breast == TRUE), TRUE, FALSE)) %>%
    rename(pval_eqtl = acat,
           pval_eqtl_breast = breast_pvalue) %>%
    select(all_of(c("group",
                    "sig_eqtl", "sig_eqtl_breast", "sig_eqtl_any", "pval_eqtl", "pval_eqtl_breast")))

  # group   pvalue  n_variants  n_features  n_indep tmi status  mtiss_acat  breast_acat n_tiss  n_features_pval_mod
  acat_sqtl <- fread(glue("../../02_acat_eqtl_sqtl/output/acat_sqtl/11tiss__{study_name}_sqtl_acat_results.tsv"))
  

  thresh_sqtl <- .05 / (acat_sqtl %>% filter(!is.na(mtiss_acat)) %>% nrow())
  thresh_sqtl_breast <- .05 / (acat_sqtl %>% filter(!is.na(breast_acat)) %>% nrow())

  acat_sqtl <- acat_sqtl %>%
    mutate(sig_sqtl = ifelse(mtiss_acat < thresh_sqtl, T, F),
           sig_sqtl_breast = ifelse(breast_acat < thresh_sqtl_breast, T, F),
           sig_sqtl_any = ifelse((sig_sqtl == TRUE | sig_sqtl_breast == TRUE), TRUE, FALSE)) %>%
    rename(pval_sqtl = mtiss_acat,
           pval_sqtl_breast = breast_acat) %>%
    # right_join(gcode_df, .,  by = c("gene_id.v26" = "group")) %>%
    select(all_of(c("group", 
                    "sig_sqtl", "sig_sqtl_breast", "sig_sqtl_any", "pval_sqtl", "pval_sqtl_breast")))

    # filter(sig_eqtl == TRUE | (sig_eqtl_breast == TRUE)) %>%
    # filter(sig_sqtl == TRUE | (sig_sqtl_breast == TRUE)) %>%
  # Have now defined our set of genes for the table, all that's left is to add/edit columns
  if (filter_sig){
    table_df <- full_join(acat_eqtl, acat_sqtl, by = c("group")) %>%
      right_join(gcode_df, .,  by = c("gene_id.v26" = "group")) %>%
      filter((sig_eqtl == TRUE | (sig_eqtl_breast == TRUE)) |
       (sig_sqtl == TRUE | (sig_sqtl_breast == TRUE)))
  } else {
    table_df <- full_join(acat_eqtl, acat_sqtl, by = c("group")) %>%
      right_join(gcode_df, .,  by = c("gene_id.v26" = "group"))
  }

  if (test){
    print(table_df %>% filter(`gene_id.v26` == 'ENSG00000198625.12') %>% pull(pval_sqtl))
    return()
  }


  # Check if we can really right join on gencode
  # if (nrow(table_df_check) != nrow(table_df)){
  #   print("Dropping rows when joining on gencode")
  # } else {
  #   print("No rows dropped with right join")
  #   return()
  # }
    
  print(glue("Table size after ACAT selection: {nrow(table_df)}"))

  # Add COJO data
  #  group  acat    breast_pvalue   group_size  group_size  available_results   available_results
  cojo_eqtl <- fread(glue("../../04_cojo_ctwas/output/condTWAS_final/{study_name}_expression_11tiss.txt")) %>%
    select(all_of(c("group", "acat", "breast_pvalue"))) %>%
    mutate(sig_eqtl_cond = ifelse(as.numeric(acat) < thresh_eqtl, T, F),
           sig_eqtl_breast_cond = ifelse(as.numeric(breast_pvalue) < thresh_eqtl_breast, T, F),
           sig_eqtl_any_cond = ifelse((sig_eqtl_cond == TRUE | sig_eqtl_breast_cond == TRUE), TRUE, FALSE)) %>%
    mutate(pval_eqtl_cond = as.numeric(acat),
           pval_eqtl_breast_cond = as.numeric(breast_pvalue)) %>%
    rename(pval_eqtl_cond_as_char = acat,
           pval_eqtl_breast_cond_as_char = breast_pvalue)


  # group   pvalue  n_variants  n_features  n_indep tmi status  mtiss_acat  breast_acat n_tiss  n_features_pval_mod cojo_status
  cojo_sqtl <- fread(glue("../../04_cojo_ctwas/output/condTWAS_final/{study_name}_11tiss.txt")) %>%
    mutate(sig_sqtl_cond = ifelse(as.numeric(mtiss_acat) < thresh_sqtl, T, F), # This columns is messy, with cojo statuses
           sig_sqtl_breast_cond = ifelse(as.numeric(breast_acat) < thresh_sqtl_breast, T, F), 
           sig_sqtl_any_cond = ifelse((sig_sqtl_cond == TRUE | sig_sqtl_breast_cond == TRUE), TRUE, FALSE)) %>%
    mutate(pval_sqtl_cond = as.numeric(mtiss_acat),
           pval_sqtl_breast_cond = as.numeric(breast_acat)) %>%
    rename(pval_sqtl_cond_as_char = mtiss_acat,
           pval_sqtl_breast_cond_as_char = breast_acat)

  table_df <- table_df %>%
    left_join(cojo_eqtl %>%
              select(all_of(c("group", "sig_eqtl_cond", "sig_eqtl_breast_cond", "sig_eqtl_any_cond",
                              "pval_eqtl_cond", "pval_eqtl_breast_cond"))),
              by = c("gene_id.v26" = "group")
              ) %>%
    left_join(cojo_sqtl %>%
              select(all_of(c("group", "sig_sqtl_cond", "sig_sqtl_breast_cond", "sig_sqtl_any_cond",
                              "pval_sqtl_cond", "pval_sqtl_breast_cond"))),
              by = c("gene_id.v26" = "group"))
  print(glue("Table size after COJO joining: {nrow(table_df)}"))

  if (study_name %in% ER_MATCH_STUDIES){
    # Add COJO data
    #  group  acat    breast_pvalue   group_size  group_size  available_results   available_results
    cojo_eqtl_er <- fread(glue("../../04_cojo_ctwas/output/condTWAS_final/{study_name}_index_ermatch_expression_11tiss.txt")) %>%
      select(all_of(c("group", "acat", "breast_pvalue"))) %>%
      mutate(sig_eqtl_cond_er = ifelse(as.numeric(acat) < thresh_eqtl, T, F),
             sig_eqtl_breast_cond_er = ifelse(as.numeric(breast_pvalue) < thresh_eqtl_breast, T, F),
             sig_eqtl_any_cond_er = ifelse((sig_eqtl_cond_er == TRUE | sig_eqtl_breast_cond_er == TRUE), TRUE, FALSE)) %>%
      mutate(pval_eqtl_cond_er = as.numeric(acat),
             pval_eqtl_breast_cond_er = as.numeric(breast_pvalue)) %>%
      rename(pval_eqtl_cond_as_char_er = acat,
             pval_eqtl_breast_cond_as_char_er = breast_pvalue)


    # group   pvalue  n_variants  n_features  n_indep tmi status  mtiss_acat  breast_acat n_tiss  n_features_pval_mod cojo_status
    cojo_sqtl_er <- fread(glue("../../04_cojo_ctwas/output/condTWAS_final/{study_name}_index_ermatch_11tiss.txt")) %>%
      mutate(sig_sqtl_cond_er = ifelse(as.numeric(mtiss_acat) < thresh_sqtl, T, F), # This columns is messy, with cojo statuses
             sig_sqtl_breast_cond_er = ifelse(as.numeric(breast_acat) < thresh_sqtl_breast, T, F), 
             sig_sqtl_any_cond_er = ifelse((sig_sqtl_cond_er == TRUE | sig_sqtl_breast_cond_er == TRUE), TRUE, FALSE)) %>%
      mutate(pval_sqtl_cond_er = as.numeric(mtiss_acat),
             pval_sqtl_breast_cond_er = as.numeric(breast_acat)) %>%
      rename(pval_sqtl_cond_as_char_er = mtiss_acat,
             pval_sqtl_breast_cond_as_char_er = breast_acat)

    table_df <- table_df %>%
      left_join(cojo_eqtl_er %>%
                select(all_of(c("group", "sig_eqtl_cond_er", "sig_eqtl_breast_cond_er", "sig_eqtl_any_cond_er",
                                "pval_eqtl_cond_er", "pval_eqtl_breast_cond_er"))),
                by = c("gene_id.v26" = "group")
                ) %>%
      left_join(cojo_sqtl_er %>%
                select(all_of(c("group", "sig_sqtl_cond_er", "sig_sqtl_breast_cond_er", "sig_sqtl_any_cond_er",
                                "pval_sqtl_cond_er", "pval_sqtl_breast_cond_er"))),
                by = c("gene_id.v26" = "group"))
    print(glue("Table size after COJO (er_match) joining: {nrow(table_df)}"))

    cojo_eqtl_ers <- fread(glue("../../04_cojo_ctwas/output/condTWAS_final/{study_name}_index_erswap_expression_11tiss.txt")) %>%
      select(all_of(c("group", "acat", "breast_pvalue"))) %>%
      mutate(sig_eqtl_cond_ers = ifelse(as.numeric(acat) < thresh_eqtl, T, F),
             sig_eqtl_breast_cond_ers = ifelse(as.numeric(breast_pvalue) < thresh_eqtl_breast, T, F),
             sig_eqtl_any_cond_ers = ifelse((sig_eqtl_cond_ers == TRUE | sig_eqtl_breast_cond_ers == TRUE), TRUE, FALSE)) %>%
      mutate(pval_eqtl_cond_ers = as.numeric(acat),
             pval_eqtl_breast_cond_ers = as.numeric(breast_pvalue)) %>%
      rename(pval_eqtl_cond_as_char_ers = acat,
             pval_eqtl_breast_cond_as_char_ers = breast_pvalue)


    # group   pvalue  n_variants  n_features  n_indep tmi status  mtiss_acat  breast_acat n_tiss  n_features_pval_mod cojo_status
    cojo_sqtl_ers <- fread(glue("../../04_cojo_ctwas/output/condTWAS_final/{study_name}_index_erswap_11tiss.txt")) %>%
      mutate(sig_sqtl_cond_ers = ifelse(as.numeric(mtiss_acat) < thresh_sqtl, T, F), # This columns is messy, with cojo statuses
             sig_sqtl_breast_cond_ers = ifelse(as.numeric(breast_acat) < thresh_sqtl_breast, T, F), 
             sig_sqtl_any_cond_ers = ifelse((sig_sqtl_cond_ers == TRUE | sig_sqtl_breast_cond_ers == TRUE), TRUE, FALSE)) %>%
      mutate(pval_sqtl_cond_ers = as.numeric(mtiss_acat),
             pval_sqtl_breast_cond_ers = as.numeric(breast_acat)) %>%
      rename(pval_sqtl_cond_as_char_ers = mtiss_acat,
             pval_sqtl_breast_cond_as_char_ers = breast_acat)

    table_df <- table_df %>%
      left_join(cojo_eqtl_ers %>%
                select(all_of(c("group", "sig_eqtl_cond_ers", "sig_eqtl_breast_cond_ers", "sig_eqtl_any_cond_ers",
                                "pval_eqtl_cond_ers", "pval_eqtl_breast_cond_ers"))),
                by = c("gene_id.v26" = "group")
                ) %>%
      left_join(cojo_sqtl_ers %>%
                select(all_of(c("group", "sig_sqtl_cond_ers", "sig_sqtl_breast_cond_ers", "sig_sqtl_any_cond_ers",
                                "pval_sqtl_cond_ers", "pval_sqtl_breast_cond_ers"))),
                by = c("gene_id.v26" = "group"))
    print(glue("Table size after COJO (er_swap) joining: {nrow(table_df)}"))
  }

  # Add enloc data (RCP)
  # group   enloc_eqtl_max_rcp_tissue   enloc_eqtl_max_locus_rcp    enloc_eqtl_lead_coloc_SNP   enloc_eqtl_lead_snp_rcpgwas_locus   locus_gwas_pip  enloc_eqtl_breast_max_locus_rcp enloc_eqtl_breast_lead_coloc_SNP    enloc_eqtl_breast_lead_snp_rcp  tissue  trait
  enloc_eqtl_fpath <- glue("../../999_reporting_friendly_tables/output/gene_rolled_enloc_eqtl_{study_name}__parse_e.txt")
  if (file.exists(enloc_eqtl_fpath)){
    have_enloc <- T
    enloc_eqtl <- fread(glue("../../999_reporting_friendly_tables/output/gene_rolled_enloc_eqtl_{study_name}__parse_e.txt")) %>%
      select(c(all_of(c("group", "enloc_eqtl_max_rcp_tissue", "enloc_eqtl_max_locus_rcp", "enloc_eqtl_min_locus_rcp",

                        "enloc_eqtl_breast_max_locus_rcp"))))

    # group   enloc_sqtl_breast_max_locus_rcp enloc_sqtl_breast_max_rcp_intron    enloc_sqtl_breast_max_lead_coloc_SNP    enloc_sqtl_breast_max_lead_snp_rcp  enloc_sqtl_max_rcp_tissue   enloc_sqtl_max_locus_rcp    enloc_sqtl_max_rcp_intron   enloc_sqtl_max_rcp_lead_coloc_SNP   enloc_sqtl_max_rcp_lead_snp_rcp intron
    enloc_sqtl <- fread(glue("../../999_reporting_friendly_tables/output/gene_rolled_enloc_sqtl_{study_name}__parse_e.txt")) %>%
      select(c(all_of(c("group", "enloc_sqtl_max_rcp_tissue", "enloc_sqtl_max_locus_rcp", "enloc_sqtl_min_locus_rcp", 
                        "enloc_sqtl_breast_max_locus_rcp"))))

    table_df <- table_df %>%
      left_join(enloc_eqtl,
                by = c("gene_id.v26" = "group")
                ) %>%
      left_join(enloc_sqtl,
                by = c("gene_id.v26" = "group"))
    print(glue("Table size after enloc joining: {nrow(table_df)}"))
  } else {
    have_enloc <- F
  }

  # Add FOCUS data (PIP)
  # ens_gene_id mol_name    ref_name    type    chrom   inference   ens_gene_id_no_decimal  num_in_focus    num_in_cred_set highest_pip highest_abs_twas_z
  focus_eqtl_fpath <- glue("../../999_reporting_friendly_tables/output/gene_rolled_focus_eqtl_11tiss_{study_name}.focus.tsv")
  if (file.exists(focus_eqtl_fpath)){
    have_focus_joint <- T
    focus_eqtl <- fread(focus_eqtl_fpath) %>%
      select(all_of(c("ens_gene_id_no_decimal", "highest_pip", "num_in_cred_set"))) %>%
      rename(gene_id_no_decimal = ens_gene_id_no_decimal,
             pip_focus_eqtl = highest_pip,
             focus_num_in_cred_set = num_in_cred_set
             )


    # ens_gene_id chrom   introns_in_focus    introns_in_cred_set highest_pip highest_abs_twas_z  max_pipmax_pip_introns  max_pip_tissues
    focus_sqtl <- fread(glue("../../999_reporting_friendly_tables/output/gene_rolled_focus_sqtl_11tiss_{study_name}.focus.tsv")) %>%
      select(all_of(c("ens_gene_id", "highest_pip", "introns_in_cred_set"))) %>%
      rename(gene_id.v26 = ens_gene_id,
             pip_focus_sqtl = highest_pip,
             focus_introns_in_cred_set = introns_in_cred_set
             )

    table_df <- table_df %>%
      left_join(focus_eqtl,
                by = c("gene_id_no_decimal" = "gene_id_no_decimal")
                ) %>%
      left_join(focus_sqtl,
                by = c("gene_id.v26" = "gene_id.v26"))

    print(glue("Table size after FOCUS joining: {nrow(table_df)}"))
  } else {
    have_focus_joint <- F
  }

  # Add FOCUS data (multi-tissue PIP columns)
  focus_eqtl_by_tiss <- fread(glue("../../999_reporting_friendly_tables/output/gene_rolled_focus_eqtl_11tiss_{study_name}_by_tissue.focus.tsv")) %>%
    select(-all_of(c("ens_gene_id", "mol_name", "type", "chrom"))) %>%
    rename_with(.cols = everything(), function(x){paste0("eqtl.", x)})

  focus_sqtl_by_tiss <- fread(glue("../../999_reporting_friendly_tables/output/gene_rolled_focus_sqtl_11tiss_{study_name}_by_tissue.focus.tsv")) %>%
    rename(`highest_pip_tissues` = `max_pip_tissues`) %>%
    rename_with(.cols = everything(), function(x){paste0("sqtl.", x)})

  table_df <- table_df %>%
    left_join(focus_eqtl_by_tiss,
              by = c("gene_id_no_decimal" = "eqtl.ens_gene_id_no_decimal")
              ) %>%
    left_join(focus_sqtl_by_tiss,
              by = c("gene_id.v26" = "sqtl.ens_gene_id"))

  print(glue("Table size after FOCUS (by tiss) joining: {nrow(table_df)}"))

  # Add CTWAS data if exists
  if (study_name %in% CTWAS_STUDIES){
    ctwas_eqtl_by_tiss <- fread(glue("../../999_reporting_friendly_tables/output/gene_rolled_ctwas_{study_name}_eqtl.susieIrss.txt")) %>%
      select(-all_of(c("chrom"))) %>%
      rename_with(.cols = everything(), function(x){paste0("ctwas.eqtl.", tolower(x))})

    ctwas_sqtl_by_tiss <- fread(glue("../../999_reporting_friendly_tables/output/gene_rolled_ctwas_{study_name}_sqtl.susieIrss.txt")) %>%
      rename_with(.cols = everything(), function(x){paste0("ctwas.sqtl.", tolower(x))})

    table_df <- table_df %>%
      left_join(ctwas_eqtl_by_tiss,
                by = c("gene_id.v26" = "ctwas.eqtl.id")
                ) %>%
      left_join(ctwas_sqtl_by_tiss,
                by = c("gene_id.v26" = "ctwas.sqtl.ens_gene_id"))

    print(glue("Table size after CTWAS (by tiss) joining: {nrow(table_df)}"))
  }

  # Focus PIP (breast only)
  focus_eqtl_breast_path <- glue("../../999_reporting_friendly_tables/output/gene_rolled_breast_only_focus_eqtl_11tiss_{study_name}.focus.tsv")
  if (file.exists(focus_eqtl_breast_path)){
    have_focus_breast <- T
    focus_eqtl_breast <- fread(focus_eqtl_breast_path) %>%
      select(all_of(c("ens_gene_id_no_decimal", "highest_pip", "num_in_cred_set"))) %>%
      rename(gene_id_no_decimal = ens_gene_id_no_decimal,
             pip_focus_breast_eqtl = highest_pip,
             focus_breast_num_in_cred_set = num_in_cred_set,
             )


    # ens_gene_id chrom   introns_in_focus    introns_in_cred_set highest_pip highest_abs_twas_z  max_pipmax_pip_introns  max_pip_tissues
    focus_sqtl_breast <- fread(glue("../../999_reporting_friendly_tables/output/gene_rolled_breast_only_focus_sqtl_11tiss_{study_name}.focus.tsv")) %>%
      select(all_of(c("ens_gene_id", "highest_pip", "introns_in_cred_set"))) %>%
      rename(gene_id.v26 = ens_gene_id,
             pip_focus_breast_sqtl = highest_pip,
             focus_breast_introns_in_cred_set = introns_in_cred_set
             )

    table_df <- table_df %>%
      left_join(focus_eqtl_breast,
                by = c("gene_id_no_decimal" = "gene_id_no_decimal")
                ) %>%
      left_join(focus_sqtl_breast,
                by = c("gene_id.v26" = "gene_id.v26"))

    print(glue("Table size after FOCUS (breast) joining: {nrow(table_df)}"))
  } else {
    have_focus_breast <- F
  }

  ## Add in cojo status
  cojo_stat_cols <- c("group", "eqtl_cojo_status", "sqtl_cojo_status", "eqtl_whitelist", "eqtl_blacklist",
                      "sqtl_whitelist", "sqtl_blacklist")
  eqtl_whitelist <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_expression_whitelist.tsv")) %>%
                      mutate(eqtl_whitelist = T,
                             eqtl_blacklist = F) %>%
                      rename(eqtl_cojo_status = status) %>%
                      select(any_of(cojo_stat_cols))
  eqtl_blacklist <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_expression_blacklist.tsv")) %>%
                      mutate(eqtl_blacklist = T) %>%
                      rename(eqtl_cojo_status = status) %>%
                      select(any_of(cojo_stat_cols))

  if (nrow(eqtl_blacklist) > 0){
    eqtl_cojo_status_df <- bind_rows(eqtl_whitelist, eqtl_blacklist) %>%
      select(any_of(cojo_stat_cols))
  } else {
    eqtl_cojo_status_df <- eqtl_whitelist
  }

  sqtl_whitelist <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_whitelist.tsv")) %>%
                      mutate(sqtl_whitelist = T,
                             sqtl_blacklist = F) %>%
                      rename(sqtl_cojo_status = status) %>%
                      select(any_of(cojo_stat_cols))
  sqtl_blacklist <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_blacklist.tsv")) %>%
                      mutate(sqtl_blacklist = T) %>%
                      rename(sqtl_cojo_status = status) %>%
                      select(any_of(cojo_stat_cols))
  if (nrow(sqtl_blacklist) > 0){
    sqtl_cojo_status_df <- bind_rows(sqtl_whitelist, sqtl_blacklist) %>%
      select(any_of(cojo_stat_cols))
  } else {
    sqtl_cojo_status_df <- sqtl_whitelist
  }

  table_df <- table_df %>%
    left_join(eqtl_cojo_status_df,
              by = c("gene_id.v26" = "group")
              ) %>%
    left_join(sqtl_cojo_status_df,
              by = c("gene_id.v26" = "group"))

  print(glue("Table size after COJO status joining: {nrow(table_df)}"))
  if (study_name %in% ER_MATCH_STUDIES){
      ## Add in cojo status
    cojo_er_stat_cols <- c("group", "eqtl_cojo_er_status", "sqtl_cojo_er_status", "eqtl_whitelist_er", "eqtl_blacklist_er",
                        "sqtl_whitelist_er", "sqtl_blacklist_er")
    eqtl_whitelist_er <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_index_ermatch_expression_whitelist.tsv")) %>%
                        mutate(eqtl_whitelist_er = T,
                               eqtl_blacklist_er = F) %>%
                        rename(eqtl_cojo_er_status = status) %>%
                        select(any_of(cojo_er_stat_cols))
    eqtl_blacklist_er <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_index_ermatch_expression_blacklist.tsv")) %>%
                        mutate(eqtl_blacklist_er = T) %>%
                        rename(eqtl_cojo_er_status = status) %>%
                        select(any_of(cojo_er_stat_cols))

    if (nrow(eqtl_blacklist_er) > 0){
      eqtl_cojo_er_status_df <- bind_rows(eqtl_whitelist_er, eqtl_blacklist_er) %>%
        select(any_of(cojo_er_stat_cols))
    } else {
      eqtl_cojo_er_status_df <- eqtl_whitelist_er
    }

    sqtl_whitelist_er <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_index_ermatch_whitelist.tsv")) %>%
                        mutate(sqtl_whitelist_er = T,
                               sqtl_blacklist_er = F) %>%
                        rename(sqtl_cojo_er_status = status) %>%
                        select(any_of(cojo_er_stat_cols))
    sqtl_blacklist_er <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_index_ermatch_blacklist.tsv")) %>%
                        mutate(sqtl_blacklist_er = T) %>%
                        rename(sqtl_cojo_er_status = status) %>%
                        select(any_of(cojo_er_stat_cols))
    if (nrow(sqtl_blacklist_er) > 0){
      sqtl_cojo_er_status_df <- bind_rows(sqtl_whitelist_er, sqtl_blacklist_er) %>%
        select(any_of(cojo_er_stat_cols))
    } else {
      sqtl_cojo_er_status_df <- sqtl_whitelist_er
    }

    table_df <- table_df %>%
      left_join(eqtl_cojo_er_status_df,
                by = c("gene_id.v26" = "group")
                ) %>%
      left_join(sqtl_cojo_er_status_df,
                by = c("gene_id.v26" = "group"))

    print(glue("Table size after COJO (ER Match) status joining: {nrow(table_df)}"))

      ## Add in cojo status
    cojo_ers_stat_cols <- c("group", "eqtl_cojo_ers_status", "sqtl_cojo_ers_status", "eqtl_whitelist_ers", "eqtl_blacklist_ers",
                        "sqtl_whitelist_ers", "sqtl_blacklist_ers")
    eqtl_whitelist_ers <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_index_erswap_expression_whitelist.tsv")) %>%
                        mutate(eqtl_whitelist_ers = T,
                               eqtl_blacklist_ers = F) %>%
                        rename(eqtl_cojo_ers_status = status) %>%
                        select(any_of(cojo_ers_stat_cols))
    eqtl_blacklist_ers <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_index_erswap_expression_blacklist.tsv")) %>%
                        mutate(eqtl_blacklist_ers = T) %>%
                        rename(eqtl_cojo_ers_status = status) %>%
                        select(any_of(cojo_ers_stat_cols))

    if (nrow(eqtl_blacklist_ers) > 0){
      eqtl_cojo_ers_status_df <- bind_rows(eqtl_whitelist_ers, eqtl_blacklist_ers) %>%
        select(any_of(cojo_ers_stat_cols))
    } else {
      eqtl_cojo_ers_status_df <- eqtl_whitelist_ers
    }

    sqtl_whitelist_ers <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_index_erswap_whitelist.tsv")) %>%
                        mutate(sqtl_whitelist_ers = T,
                               sqtl_blacklist_ers = F) %>%
                        rename(sqtl_cojo_ers_status = status) %>%
                        select(any_of(cojo_ers_stat_cols))
    sqtl_blacklist_ers <- fread(glue("../../04_cojo_ctwas/output/intermediate_data/cojo_input/{study_name}_index_erswap_blacklist.tsv")) %>%
                        mutate(sqtl_blacklist_ers = T) %>%
                        rename(sqtl_cojo_ers_status = status) %>%
                        select(any_of(cojo_ers_stat_cols))
    if (nrow(sqtl_blacklist_ers) > 0){
      sqtl_cojo_ers_status_df <- bind_rows(sqtl_whitelist_ers, sqtl_blacklist_ers) %>%
        select(any_of(cojo_ers_stat_cols))
    } else {
      sqtl_cojo_ers_status_df <- sqtl_whitelist_ers
    }

    table_df <- table_df %>%
      left_join(eqtl_cojo_ers_status_df,
                by = c("gene_id.v26" = "group")
                ) %>%
      left_join(sqtl_cojo_ers_status_df,
                by = c("gene_id.v26" = "group"))

    print(glue("Table size after COJO (ER Swap) status joining: {nrow(table_df)}"))
  }


	# COLOC
	coloc_eqtl_path <- glue("../../09_coloc/output/final/gene_rolled_eqtl_{tiss_list_basename_no_ext}_{study_name}.tsv")
	coloc_sqtl_path <- glue("../../09_coloc/output/final/gene_rolled_sqtl_{tiss_list_basename_no_ext}_{study_name}.tsv")
	if (file.exists(glue(coloc_eqtl_path))){
	  coloc_eqtl_df <- fread(coloc_eqtl_path) %>%
      rename_with(.cols = everything(), function(x){paste0("eqtl.", x)})
	  coloc_sqtl_df <- fread(coloc_sqtl_path) %>%
      rename_with(.cols = everything(), function(x){paste0("sqtl.", x)})

		table_df <- table_df %>%
			left_join(coloc_eqtl_df, by = c("gene_id.v26" = "eqtl.gene")) %>%
			left_join(coloc_sqtl_df, by = c("gene_id.v26" = "sqtl.gene"))
	}

	# Check if coloc names exist
	coloc_names <- c()
	if ("eqtl.coloc_top_h4_tissue" %in% names(table_df)){
	  coloc_names <- c(coloc_names, "eqtl.coloc_top_h4", "eqtl.coloc_top_h4_tissue", "eqtl.coloc_top_h4_snp", "eqtl.coloc_top_h4_snp_post")
	}
	if ("sqtl.coloc_top_h4_tissue" %in% names(table_df)){
		coloc_names <- c(coloc_names, "sqtl.coloc_top_h4", "sqtl.coloc_top_h4_tissue",
										 "sqtl.coloc_top_h4_snp", "sqtl.coloc_top_h4_snp_post", "sqtl.coloc_top_h4_intron")
	}

  ## Add in zscore columns
  table_df <- table_df %>% add_zscore_tissue(study_name,
                                             unique(table_df$gene_id.v26),
                                             "eqtl" , tiss_filter)
  table_df <- table_df %>% add_zscore_tissue(study_name,
                                             unique(table_df$gene_id.v26),
                                             "sqtl" , tiss_filter)
  print(glue("Table size after adding zscore columns: {nrow(table_df)}"))

  ## Can return early for use with Guimin wanted tables
  table_df <- table_df %>%
    mutate(thresh_eqtl = thresh_eqtl,
           thresh_sqtl = thresh_sqtl,
           thresh_eqtl_breast = thresh_eqtl_breast,
           thresh_sqtl_breast = thresh_sqtl_breast)
  if (return_raw == T){
    print("Returning raw combined results early.")
    return(table_df)
  }

  ## Create columns for display
  final_col_order <- c("Gene Loci", "gene_name.v40", "gene_type.v40", "gene_type.v26", "gene_name.v38", "gene_id_no_decimal", "gene_id.v26", "genomic_coord", "genomic_coord.v26", "Index SNP", "Distance from Index SNP(kb)",
                      
                       # These cols expand to eqtl/sqtl specific columns
                       "min_p_11_tiss",
                       "pval_eqtl", "pval_sqtl",
                       
                       "min_p_breast",
                       "pval_eqtl_breast", "pval_sqtl_breast",

                       "min_cojo_p_11_tiss",
                       "min_cojo_er_p_11_tiss",
                       "min_cojo_ers_p_11_tiss",
                       "pval_eqtl_cond", "pval_sqtl_cond",
                       "pval_eqtl_cond_er", "pval_sqtl_cond_er",
                       "pval_eqtl_cond_ers", "pval_sqtl_cond_ers",

                       "min_cojo_p_breast",
                       "min_cojo_er_p_breast",
                       "min_cojo_ers_p_breast",
                       "pval_eqtl_breast_cond", "pval_sqtl_breast_cond",
                       "pval_eqtl_breast_cond_er", "pval_sqtl_breast_cond_er",
                       "pval_eqtl_breast_cond_ers", "pval_sqtl_breast_cond_ers",
                       ## END expanded columns
                       
                       "max_rcp_11_tiss", "max_rcp_tissue",
                       "max_rcp_breast",

                       "max_pip_11_tiss", "max_pip_breast",

                       "twas_approach_id",
                       "ctwas_approach_id",
                       "ctwas_er_approach_id",
                       "ctwas_ers_approach_id",

                       # Asked for 1/23/22
                       "eqtl_max_abs_zscore_tiss",
                       "sqtl_max_abs_zscore_tiss", 

                       # FOCUS by tissue target PIP columns
                      "adipose_subcutaneous_max_pip",
                      "adipose_subcutaneous_max_susie_pip",
                      "adipose_visceral_omentum_max_pip",
                      "adipose_visceral_omentum_max_susie_pip",
                      "breast_mammary_tissue_max_pip",
                      "breast_mammary_tissue_max_susie_pip",
                      "cells_cultured_fibroblasts_max_pip",
                      "cells_cultured_fibroblasts_max_susie_pip",
                      "cells_ebv-transformed_lymphocytes_max_pip",
                      "cells_ebv-transformed_lymphocytes_max_susie_pip",
                      "liver_max_pip",
                      "liver_max_susie_pip",
                      "ovary_max_pip",
                      "ovary_max_susie_pip",
                      "spleen_max_pip",
                      "spleen_max_susie_pip",
                      "uterus_max_pip",
                      "uterus_max_susie_pip",
                      "vagina_max_pip",
                      "vagina_max_susie_pip",
                      "whole_blood_max_pip",
                      "whole_blood_max_susie_pip",

                      # max eqtl pip/tissues
                      # max sqtl pip/tissues
                      "eqtl.highest_pip",
                      "ctwas.eqtl.highest_susie_pip",
                      "eqtl.highest_pip_tissues",
                      "ctwas.eqtl.highest_susie_pip_tissues",
                      "sqtl.highest_pip",
                      "ctwas.sqtl.highest_susie_pip",
                      "sqtl.highest_pip_tissues",
                      "ctwas.sqtl.highest_susie_pip_tissues",

											coloc_names,

                       # Not asked for but potentially useful
                       "gene_id.v26",
                       "thresh_eqtl", "thresh_sqtl",
                       "thresh_eqtl_breast", "thresh_sqtl_breast",

                       "eqtl_cojo_status", "sqtl_cojo_status",
                       "eqtl_cojo_er_status", "sqtl_cojo_er_status",
                       "eqtl_cojo_ers_status", "sqtl_cojo_ers_status",
                       "eqtl_whitelist", "eqtl_blacklist",
                       "sqtl_whitelist", "sqtl_blacklist",
                       "eqtl_whitelist_er", "eqtl_blacklist_er",
                       "sqtl_whitelist_er", "sqtl_blacklist_er",
                       "eqtl_whitelist_ers", "eqtl_blacklist_ers",
                       "sqtl_whitelist_ers", "sqtl_blacklist_ers"
                       )
  table_df <- table_df %>%
    mutate(min_p_11_tiss = pmin(pval_eqtl, pval_sqtl, na.rm=T),
           min_p_breast = pmin(pval_eqtl_breast, pval_sqtl_breast, na.rm=T),
           min_cojo_p_11_tiss = pmin(pval_eqtl_cond, pval_sqtl_cond, na.rm=T),
           min_cojo_p_breast = pmin(pval_eqtl_breast_cond, pval_sqtl_breast_cond, na.rm=T)
           ) %>%
    mutate(twas_approach_eqtl = ifelse(sig_eqtl_any == TRUE, "express", ""),
           twas_approach_sqtl = ifelse(sig_sqtl_any == TRUE, "splice", ""),
           twas_approach_id = glue("{twas_approach_eqtl}|{twas_approach_sqtl}")) %>%
    mutate(twas_approach_id = str_replace(twas_approach_id, "NA", ""),
           twas_approach_id = str_replace(twas_approach_id, "^\\|", ""),
           twas_approach_id = str_replace(twas_approach_id, "\\|$", "")) %>%
    mutate(ctwas_approach_eqtl = ifelse(sig_eqtl_any_cond == TRUE, "express", ""),
           ctwas_approach_sqtl = ifelse(sig_sqtl_any_cond == TRUE, "splice", ""),
           ctwas_approach_id = glue("{ctwas_approach_eqtl}|{ctwas_approach_sqtl}")) %>%
    mutate(ctwas_approach_id = str_replace(ctwas_approach_id, "NA", ""),
           ctwas_approach_id = str_replace(ctwas_approach_id, "^\\|", ""),
           ctwas_approach_id = str_replace(ctwas_approach_id, "\\|$", ""),
           ctwas_approach_id = str_replace(ctwas_approach_id, "NA", "")
           )
  if (have_enloc){
    table_df <- table_df %>%
      mutate(
           max_rcp_11_tiss = pmax(enloc_eqtl_max_locus_rcp, enloc_sqtl_max_locus_rcp, na.rm=T),
           max_rcp_tissue = ifelse(enloc_eqtl_max_locus_rcp == max_rcp_11_tiss, enloc_eqtl_max_rcp_tissue, 
                                   enloc_sqtl_max_rcp_tissue),
           max_rcp_breast = pmax(enloc_eqtl_breast_max_locus_rcp, enloc_sqtl_breast_max_locus_rcp, na.rm=T))
  }

  if (have_focus_joint){
    table_df <- table_df %>%
      mutate(
           max_pip_11_tiss = pmax(pip_focus_eqtl, pip_focus_sqtl, na.rm=T))
  }
  if (have_focus_breast){
    table_df <- table_df %>%
      mutate(
           max_pip_breast = pmax(pip_focus_breast_eqtl, pip_focus_breast_sqtl, na.rm=T))
  }


  if (study_name %in% ER_MATCH_STUDIES){
    table_df <- table_df %>%
      mutate(
           min_cojo_er_p_11_tiss = pmin(pval_eqtl_cond_er, pval_sqtl_cond_er, na.rm=T),
           min_cojo_er_p_breast = pmin(pval_eqtl_breast_cond_er, pval_sqtl_breast_cond_er, na.rm=T)
      ) %>%
      mutate(ctwas_er_approach_eqtl = ifelse(sig_eqtl_any_cond_er == TRUE, "express", ""),
             ctwas_er_approach_sqtl = ifelse(sig_sqtl_any_cond_er == TRUE, "splice", ""),
             ctwas_er_approach_id = glue("{ctwas_er_approach_eqtl}|{ctwas_er_approach_sqtl}")) %>%
      mutate(ctwas_er_approach_id = str_replace(ctwas_er_approach_id, "NA", ""),
             ctwas_er_approach_id = str_replace(ctwas_er_approach_id, "^\\|", ""),
             ctwas_er_approach_id = str_replace(ctwas_er_approach_id, "\\|$", ""),
             ctwas_er_approach_id = str_replace(ctwas_er_approach_id, "NA", "")
             )

    table_df <- table_df %>%
      mutate(
           min_cojo_ers_p_11_tiss = pmin(pval_eqtl_cond_ers, pval_sqtl_cond_ers, na.rm=T),
           min_cojo_ers_p_breast = pmin(pval_eqtl_breast_cond_ers, pval_sqtl_breast_cond_ers, na.rm=T)) %>%
      mutate(ctwas_ers_approach_eqtl = ifelse(sig_eqtl_any_cond_ers == TRUE, "express", ""),
             ctwas_ers_approach_sqtl = ifelse(sig_sqtl_any_cond_ers == TRUE, "splice", ""),
             ctwas_ers_approach_id = glue("{ctwas_ers_approach_eqtl}|{ctwas_ers_approach_sqtl}")) %>%
      mutate(ctwas_ers_approach_id = str_replace(ctwas_ers_approach_id, "NA", ""),
             ctwas_ers_approach_id = str_replace(ctwas_ers_approach_id, "^\\|", ""),
             ctwas_ers_approach_id = str_replace(ctwas_ers_approach_id, "\\|$", ""),
             ctwas_ers_approach_id = str_replace(ctwas_ers_approach_id, "NA", "")
             )
  }

  table_df <- table_df %>%
    mutate(
      adipose_subcutaneous_max_pip = pmax(eqtl.adipose_subcutaneous_highest_pip,
                                               sqtl.adipose_subcutaneous_highest_pip, na.rm=T),
      adipose_visceral_omentum_max_pip = pmax(eqtl.adipose_visceral_omentum_highest_pip,
                                               sqtl.adipose_visceral_omentum_highest_pip, na.rm=T),
      breast_mammary_tissue_max_pip = pmax(eqtl.breast_mammary_tissue_highest_pip,
                                               sqtl.breast_mammary_tissue_highest_pip, na.rm=T),
      cells_cultured_fibroblasts_max_pip = pmax(eqtl.cells_cultured_fibroblasts_highest_pip,
                                               sqtl.cells_cultured_fibroblasts_highest_pip, na.rm=T),
      `cells_ebv-transformed_lymphocytes_max_pip` = pmax(`eqtl.cells_ebv-transformed_lymphocytes_highest_pip`,
                                               `sqtl.cells_ebv-transformed_lymphocytes_highest_pip`, na.rm=T),
      liver_max_pip = pmax(eqtl.liver_highest_pip,
                                               sqtl.liver_highest_pip, na.rm=T),
      ovary_max_pip = pmax(eqtl.ovary_highest_pip,
                                               sqtl.ovary_highest_pip, na.rm=T),
      spleen_max_pip = pmax(eqtl.spleen_highest_pip,
                                               sqtl.spleen_highest_pip, na.rm=T),
      uterus_max_pip = pmax(eqtl.uterus_highest_pip,
                                               sqtl.uterus_highest_pip, na.rm=T),
      vagina_max_pip = pmax(eqtl.vagina_highest_pip,
                                               sqtl.vagina_highest_pip, na.rm=T),
      whole_blood_max_pip = pmax(eqtl.whole_blood_highest_pip,
                                               sqtl.whole_blood_highest_pip, na.rm=T))
  if (study_name %in% CTWAS_STUDIES){
    table_df <- table_df %>%
      mutate(
        adipose_subcutaneous_max_susie_pip = pmax(ctwas.eqtl.adipose_subcutaneous_highest_susie_pip,
                                                 ctwas.sqtl.adipose_subcutaneous_highest_susie_pip, na.rm=T),
        adipose_visceral_omentum_max_susie_pip = pmax(ctwas.eqtl.adipose_visceral_omentum_highest_susie_pip,
                                                 ctwas.sqtl.adipose_visceral_omentum_highest_susie_pip, na.rm=T),
        breast_mammary_tissue_max_susie_pip = pmax(ctwas.eqtl.breast_mammary_tissue_highest_susie_pip,
                                                 ctwas.sqtl.breast_mammary_tissue_highest_susie_pip, na.rm=T),
        cells_cultured_fibroblasts_max_susie_pip = pmax(ctwas.eqtl.cells_cultured_fibroblasts_highest_susie_pip,
                                                 ctwas.sqtl.cells_cultured_fibroblasts_highest_susie_pip, na.rm=T),
        `cells_ebv-transformed_lymphocytes_max_susie_pip` = pmax(`ctwas.eqtl.cells_ebv-transformed_lymphocytes_highest_susie_pip`,
                                                 `ctwas.sqtl.cells_ebv-transformed_lymphocytes_highest_susie_pip`, na.rm=T),
        liver_max_susie_pip = pmax(ctwas.eqtl.liver_highest_susie_pip,
                                                 ctwas.sqtl.liver_highest_susie_pip, na.rm=T),
        ovary_max_susie_pip = pmax(ctwas.eqtl.ovary_highest_susie_pip,
                                                 ctwas.sqtl.ovary_highest_susie_pip, na.rm=T),
        spleen_max_susie_pip = pmax(ctwas.eqtl.spleen_highest_susie_pip,
                                                 ctwas.sqtl.spleen_highest_susie_pip, na.rm=T),
        uterus_max_susie_pip = pmax(ctwas.eqtl.uterus_highest_susie_pip,
                                                 ctwas.sqtl.uterus_highest_susie_pip, na.rm=T),
        vagina_max_susie_pip = pmax(ctwas.eqtl.vagina_highest_susie_pip,
                                                 ctwas.sqtl.vagina_highest_susie_pip, na.rm=T),
        whole_blood_max_susie_pip = pmax(ctwas.eqtl.whole_blood_highest_susie_pip,
                                                 ctwas.sqtl.whole_blood_highest_susie_pip, na.rm=T),
        max_susie_pip_11_tiss = pmax(ctwas.eqtl.highest_susie_pip, ctwas.sqtl.highest_susie_pip, na.rm=T)
      )
  }
  if (!is.na(table_names[1])){
    write_delim(table_df %>% select(any_of(final_col_order)), file.path(odir, glue(table_names[1])), delim="\t")
  }

  # # Filter for different table subtypes

  # ENSG	gene_name38	gene_type38	chrnum	strand38	pos38_start	pos38_end	hg38_position	pheno_overall	pheno_ERpos	pheno_ERneg	FirstAuthors	PUBMEDID	subsource	publication_date_first	pvalue_overall	er_pos_pvalue	er_neg_pvalue	Z_ov_Kar	P_ov_Kar	Z_ERn_Kar	P_ERn_Kar	Z_ov_Wen	P_ov_Wen	Z_ov_Jia	P_ov_Jia	Z_ERp_Jia	P_ERp_Jia	Z_ERn_Jia	P_ERn_Jia	note
  prev_twas_df <- fread("../input/Ten_TWAS_papers_gene_list_2023Feb28.tsv")

  # Genes not reported in ../input/Ten_TWAS_papers_gene_list_2023Feb28.tsv
  table_df_twas_new <- table_df %>%
    filter(!(gene_id_no_decimal %in% prev_twas_df$ENSG))
  print(glue("{nrow(table_df_twas_new)} genes out of {nrow(table_df)} not previously reported by TWAS"))
  if (!is.na(table_names[2])){
    write_delim(table_df_twas_new %>% select(any_of(final_col_order)), file.path(odir, glue(table_names[2])), delim="\t")
  }

  # Still significant (at the same threshold) after COJO conditioning
  if (study_name %in% ER_MATCH_STUDIES){
    cojo_sig_check_df <- table_df %>%
      filter((sig_eqtl_any_cond == TRUE) | (sig_sqtl_any_cond == TRUE))
    cojo_sig_df <- table_df %>%
      filter((sig_eqtl_any_cond == TRUE) | (sig_sqtl_any_cond == TRUE) |
             (sig_eqtl_any_cond_er == TRUE) | (sig_sqtl_any_cond_er == TRUE) | (sig_eqtl_any_cond_ers) | (sig_sqtl_any_cond_ers == TRUE))
    print(glue("{nrow(cojo_sig_df)} genes out of {nrow(table_df)} still significant after COJO (er match)"))
    print(glue("{nrow(cojo_sig_check_df)} genes out of {nrow(table_df)} would still be significant after non-er match/swap COJO only."))
  } else {
    cojo_sig_df <- table_df %>%
      filter((sig_eqtl_any_cond == TRUE) | (sig_sqtl_any_cond == TRUE))
    print(glue("{nrow(cojo_sig_df)} genes out of {nrow(table_df)} still significant after COJO"))
  }
  if (!is.na(table_names[3])){
    write_delim(cojo_sig_df %>% select(any_of(final_col_order)), file.path(odir, glue(table_names[3])), delim="\t")
  }

  if (!is.na(table_names[4])){
    cand_caus_df <- filter_for_caus_cand(study_name, table_df, gcode_df, "gene_id.v26", "gene_id_no_decimal")
    write_delim(cand_caus_df %>% select(any_of(final_col_order)), file.path(odir, glue(table_names[4])), delim="\t")
    # In at least 1 FOCUS credible set. For tissues in cred set, RCP must be >
    # 0.5 for all intron/gene records in all tissues for enloc.
    print(glue("{nrow(cand_caus_df)} genes out of {nrow(table_df)} are in at least 1 FOCUS cred set and have RCP > 0.5 for all records in those cred set tissues."))
  }


  if (!is.na(table_names[5])){
    ## Additional tables
    ## Supp Table 2/3
    # Everything in table 1 ER- Gene, Tissue, ZScore
    new_table_cols <- c("gene_name.v40", "gene_type.v40", "tissue", "zscore", "gene_id.v26", "gene_type.v26")

    eqtl_tiss_gene_z_df <- left_join(table_df %>% 
                                     select(any_of(new_table_cols)),
                                     get_zscore_detail(study_name, table_df$gene_id.v26, "eqtl", "table", tiss_filter),
                                     by= c("gene_id.v26")) %>%
      select(all_of(new_table_cols)) %>%
      arrange(gene_name.v40, tissue)
      
    sqtl_tiss_gene_z_df <- left_join(table_df %>% 
                                     select(any_of(new_table_cols)),
                                     get_zscore_detail(study_name, table_df$gene_id.v26, "sqtl", "table", tiss_filter),
                                     by= c("gene_id.v26")) %>%
      select(all_of(new_table_cols)) %>%
      arrange(gene_name.v40, tissue)
    eqtl_name_pat <- glue("eqtl_{glue('{table_names[5]}')}")
    sqtl_name_pat <- glue("sqtl_{glue('{table_names[5]}')}")

    write_delim(eqtl_tiss_gene_z_df, file.path(odir, glue(eqtl_name_pat)), delim="\t")
    write_delim(sqtl_tiss_gene_z_df, file.path(odir, glue(sqtl_name_pat)), delim="\t")
  }

  # Everything in table 1 ER+ Gene, Tissue, ZScore

  print(glue(""))
  return()
}
# In at least 1 FOCUS credible set. For tissues in cred set, RCP must be > 0.5
# for all intron/gene records in all tissues for enloc.
filter_for_caus_cand <- function(study_name, in_df, gcode_df, ensg_v26_full_col, ensg_v26_no_dec_col){
  require(tidyverse)
  require(data.table)
  require(glue)
  require(purrr)

  GROUPINGS <- fread("../input/combine_phenotype_groups.txt.gz") %>% 
    as_tibble() %>% 
    rename(group = gene_id) %>%
    filter(group %in% in_df[[ensg_v26_full_col]])
  GROUPINGS[c("tissue", "intron")] <- stringr::str_split_fixed(GROUPINGS$intron_id, pattern="\\.", 2)
  GROUPINGS$tissue_og <- GROUPINGS$tissue
  GROUPINGS$tissue <- tolower(GROUPINGS$tissue)

  # Read in raw FOCUS data
  focus_want <- c("gene_id.v26", "tissue")
  focus_eqtl <- fread(glue("../../05_focus/output/focus_eqtl_11tiss_{study_name}.focus.tsv")) %>%
    filter(ens_gene_id != "NULL.MODEL") %>% 
    mutate(ens_gene_id_no_decimal = str_extract(ens_gene_id, "(ENSG\\d{1,})")) %>%
    filter(ens_gene_id_no_decimal %in% in_df[[ensg_v26_no_dec_col]],
           in_cred_set == 1) %>%
    inner_join(gcode_df %>% 
               select(all_of(c("gene_id.v26", "gene_id_no_decimal"))),
               by = c("ens_gene_id_no_decimal" = "gene_id_no_decimal")) %>%
    select(any_of(focus_want))

  focus_sqtl <- fread(glue("../../05_focus/output/focus_sqtl_11tiss_{study_name}.focus.tsv")) %>%
    filter(ens_gene_id != "NULL.MODEL",
           in_cred_set == 1) %>%
    rename(intron = ens_gene_id) %>%
    inner_join(GROUPINGS, by = c("intron", "tissue")) %>%
    rename(!!ensg_v26_full_col := group) %>%
    select(-c("intron_id", "gtex_intron_id")) %>%
    select(any_of(focus_want)) %>%
    distinct()

  focus_cred_set_tissues <- bind_rows(focus_eqtl, focus_sqtl) %>% distinct() %>% arrange(gene_id.v26)
  focus_cred_set_ensg <- focus_cred_set_tissues %>% pull(ensg_v26_full_col) %>% unique()

  # Read in raw enloc data
  enloc_eqtl <- fread(glue("../../03_enloc/output/final_enloc_result_eqtl/{study_name}__parse_e.txt")) %>%
    filter(molecular_qtl_trait %in% in_df[[ensg_v26_full_col]]) %>%
    rename(!!ensg_v26_full_col := molecular_qtl_trait)
  enloc_eqtl$tissue <- tolower(enloc_eqtl$tissue)

  enloc_sqtl <- fread(glue("../../03_enloc/output/final_enloc_result_sqtl/{study_name}__parse_e.txt")) %>%
    rename(intron = molecular_qtl_trait)
  enloc_sqtl$tissue <- tolower(enloc_sqtl$tissue)
  enloc_sqtl <- enloc_sqtl %>%
    inner_join(., GROUPINGS, by = c("intron", "tissue")) %>%
    select(-intron_id) %>%
    rename(!!ensg_v26_full_col := group)


  check_if_gene_caus_cand <- function(cur_gene, focus_cred_set_tissues, enloc_eqtl, enloc_sqtl){
    require(tidyverse)
    # if (cur_gene == "ENSG00000198807.12"){
    #   browser()
    # }
    # Get the tissues the cur_gene is in the cred set for
    check_tiss <- focus_cred_set_tissues %>%
      filter_at(ensg_v26_full_col, all_vars(. == cur_gene)) %>%
      pull(tissue)
    tiss_check_status <- c() # T for pass, F for not
    for (cur_tiss in check_tiss){
      min_eqtl_locus_rcp <- enloc_eqtl %>%
        filter_at(ensg_v26_full_col, all_vars(. == cur_gene)) %>%
        filter(tissue == cur_tiss) %>%
        pull(locus_rcp) %>% min(na.rm=T)

      min_sqtl_locus_rcp <- enloc_sqtl %>%
        filter_at(ensg_v26_full_col, all_vars(. == cur_gene)) %>%
        filter(tissue == cur_tiss) %>%
        pull(locus_rcp) %>% min(na.rm=T)

      tiss_check_status <- c(tiss_check_status, all(c(min_eqtl_locus_rcp, min_sqtl_locus_rcp) > 0.5))
    }
    return(tibble(rcp_pass=all(tiss_check_status)))
  
  }
  focus_cred_set_ensg_rcp_checked <- bind_cols(tibble(gene_id.v26 = focus_cred_set_ensg),
                                               focus_cred_set_ensg %>%
                                                 map_dfr(~check_if_gene_caus_cand(.,
                                                         focus_cred_set_tissues,
                                                         enloc_eqtl, enloc_sqtl)))
  out_df <- left_join(in_df, focus_cred_set_ensg_rcp_checked) %>%
    filter(rcp_pass == TRUE)
    
}

add_zscore_tissue <- function(base_table, study_name, ensg_v26_genes, eqtl_sqtl, tiss_filter){
  zscore_tissue <- get_zscore_detail(study_name, ensg_v26_genes, eqtl_sqtl, want="tissue", tiss_filter) %>%
    select(-all_of(c("zscore_for_max_z", "pval_for_max_z"))) %>% distinct()
    
  rdf <- left_join(base_table, zscore_tissue, by = c("gene_id.v26"))
  return(rdf)
}

# `study_name` - str, study ID
# `engs_v26_genes` - vector of gencode 26 ENSG with decimal genes to filter for
# `eqtl_sqtl` - Investigating "eqtl" or "sqtl" based zscores?
# `want` - Want a "table" or simply the "tissue" with largest absolute z-score?
#        `table` - Returns tibble with cols: `gene_id.v26`, `tissue`, `zscore`
#        `tissue` 
get_zscore_detail <- function(study_name, ensg_v26_genes, eqtl_sqtl, want, tiss_filter){
  require(tidyverse)
  require(data.table)
  require(glue)
  require(purrr)

  table_df <- get_zscore_table(study_name, ensg_v26_genes, eqtl_sqtl, tiss_filter)
  if (want == "table"){
    return(table_df)
  } else if (want == "tissue") {
    max_abs_z_df <- table_df %>%
      mutate(abs_zscore = abs(zscore)) %>%
      group_by(gene_id.v26) %>%
      summarize(max_abs_zscore = max(abs_zscore, na.rm=T))
    
    get_one_gene_max_abs_z_tiss <- function(cur_gene){
      cur_gene_df <- table_df %>%
        filter(gene_id.v26 == cur_gene) %>%
        mutate(abs_zscore = abs(zscore))
      cur_gene_abs_key <- max_abs_z_df %>%
        filter(gene_id.v26 == cur_gene)

      cur_max_abs_z_tiss <- cur_gene_df %>%
        filter(abs_zscore == cur_gene_abs_key$max_abs_zscore)

      pval_for_max_z <- cur_max_abs_z_tiss %>% pull(pvalue) %>% unique()
      zscore_for_max_z <- cur_max_abs_z_tiss %>% pull(zscore) %>% unique()
      abs_zscore_for_max_z <- cur_max_abs_z_tiss %>% pull(abs_zscore) %>% unique()

      if (length(pval_for_max_z) > 1){
        if (length(abs_zscore_for_max_z) == 1){
          pval_for_max_z <- pval_for_max_z

        }
      }
      if (length(zscore_for_max_z) > 2){
        browser()
      }
      if (length(zscore_for_max_z) > 1){
        zscore_for_max_z <- glue("+-{abs_zscore_for_max_z}")
      }

      if (nrow(cur_max_abs_z_tiss) > 1){
        combined_tiss_names <- paste(unique(cur_max_abs_z_tiss$tissue), collapse=", ")
        rtib <- tibble(gene_id.v26 = cur_gene,
                       !!glue("{eqtl_sqtl}_max_abs_zscore_tiss") := combined_tiss_names,
                       zscore_for_max_z = zscore_for_max_z %>% as.character(),
                       pval_for_max_z = pval_for_max_z)
      } else  {
        rtib <- tibble(gene_id.v26 = cur_gene,
                       !!glue("{eqtl_sqtl}_max_abs_zscore_tiss") := cur_max_abs_z_tiss$tissue,
                       zscore_for_max_z = zscore_for_max_z %>% as.character(),
                       pval_for_max_z = pval_for_max_z)
      }
      return(rtib)
    }
    return_tissue_df <- ensg_v26_genes %>% map_dfr(~get_one_gene_max_abs_z_tiss(.))
    return(return_tissue_df)
  }
}

get_zscore_table <- function(study_name, ensg_v26_genes, eqtl_sqtl, tiss_filter){
  require(tidyverse)
  require(data.table)
  require(glue)
  final_table_cols <- c("gene_id.v26", "tissue", "zscore", "pvalue")
  if (eqtl_sqtl == "eqtl"){
    by_tiss_files <- list.files(path="../../01_spredixcan_eqtl_sqtl/output/spredixcan_eqtl_mashr/", pattern="spredixcan_igwas_gtexmashrv8_meta_analysis_BCAC_UKB_ovr_brca*", full.names=T)
    # gene,gene_name,zscore,effect_size,pvalue,var_g,pred_perf_r2,pred_perf_pval,pred_perf_qval,n_snps_used,n_snps_in_cov,n_snps_in_model,best_gwas_p,largest_weight

    master <- tibble()
    if (!is.na(tiss_filter)){
      tiss_vec <- fread(tiss_filter, header=F) %>% pull(V1)
    }
    for (cur_file in by_tiss_files){
      cur_file_base <- basename(cur_file)
      tissue_name <- str_match(cur_file_base, "__PM__(.*).csv")[2]

      if (!is.na(tiss_filter)){
        if (!(tissue_name %in% tiss_vec)){
          next
        }
      }
      cur_df <- fread(cur_file)  %>%
        filter(gene %in% ensg_v26_genes) %>%
        mutate(tissue = tissue_name) %>%
        rename(gene_id.v26 = gene) %>%
        select(all_of(final_table_cols))
      
      master <- bind_rows(master, cur_df)
    }
    master <- master %>%
      arrange(gene_id.v26, tissue)
  } else if (eqtl_sqtl == "sqtl"){
    # Key for introns to genes
    GROUPINGS <- fread("../input/combine_phenotype_groups.txt.gz") %>% 
      as_tibble() %>% 
      rename(group = gene_id) %>%
      filter(group %in% ensg_v26_genes)
    GROUPINGS[c("tissue", "intron")] <- stringr::str_split_fixed(GROUPINGS$intron_id, pattern="\\.", 2)

    master <- fread(glue("../../02_acat_eqtl_sqtl/input/acat_prep_sqtl/{study_name}_sqtl_acat_input.csv")) %>%
      as_tibble()

    master[c("tissue", "intron")] <- stringr::str_split_fixed(master$gene, pattern="\\.", 2)
    master <- master %>% inner_join(GROUPINGS, by = c("intron", "tissue")) %>%
      rename(gene_id.v26 = group) %>%
      select(all_of(final_table_cols)) %>%
      arrange(gene_id.v26, tissue)
    if (!is.na(tiss_filter)){
      tiss_vec <- fread(tiss_filter, header=F) %>% pull(V1)
      master <- master %>%
        filter(tissue %in% tiss_vec)
    }
  } else {
    stop("eqtl_sqtl needs to be 'eqtl' or 'sqtl'")
  }
  return(master)
}


if (!interactive()) {
  # STUDIES %>% walk(~make_tables(.))
  # i_subtype_studies <- paste0("intrinsic_subtype_", 1:5)
  # i_subtype_studies %>% walk(~make_tables(., table_names=c("supplementary_table_{study_name}.tsv"), filter_sig=F))
  # make_tables("meta_analysis_BCAC_UKB_ovr_brca")
} else {
  # STUDIES %>% walk(~make_tables(.))
  # i_subtype_studies <- paste0("intrinsic_subtype_", 1:5)
  # i_subtype_studies %>% walk(~make_tables(., table_names=c("supplementary_table_{study_name}.tsv"), filter_sig=F))
  # make_tables("meta_analysis_BCAC_UKB_ovr_brca")
}
