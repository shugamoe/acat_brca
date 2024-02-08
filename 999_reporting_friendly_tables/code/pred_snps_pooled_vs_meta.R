
main <- function(
  study_name = "meta_analysis_BCAC_UKB_ovr_brca",
  tiss_filter = "../../02_acat_eqtl_sqtl/input/tissue_lists/11tiss.txt",
  gcode_convert = "../output/gencode_26_40_conversion_and_distances_to_known_brca_index_snps.tsv",
  odir = "../output/pooled_vs_meta/"
){
  require(tidyverse)
  require(data.table)
  require(glue)
  source("dz_jl_tables_for_jcm.R")

  table_df <- runonce::save_run(make_tables(study_name=study_name,
                         tiss_filter=tiss_filter,
                         gcode_convert=gcode_convert,
                         return_raw=T), glue("../output/tmp/{study_name}_table.rds"))

  table_filt_df <- table_df %>%
    filter(sig_eqtl_any == T | sig_sqtl_any == T)
  names_filt_df <- table_filt_df %>%
    select(all_of(c("gene_id.v26", "gene_name.v26", "gene_name.v40")))

  # Retrieve meta ukb
  meta_ukb_df <- fread("../input/Anno_res_metaBC.allChr.breast.logistic.gz") %>%
    select(all_of(c("rsID", "beta_meta", "SE_meta", "p_meta", "REF", "ALT"))) %>%
    rename(REF_meta = REF,
           ALT_meta = ALT)

  # Retrieve inc vs prev pooled ukb
  pooled_ukb_df <- fread("../input/Anno_res_anyBC_full_chrom.breast.logistic.gz") %>%
    select(all_of(c("rsID", "beta_ALT", "SE", "P", "REF", "ALT"))) %>%
    rename(beta_pool = beta_ALT,
           SE_pool = SE,
           p_pool = P,
           REF_pool = REF,
           ALT_pool = ALT
           )

  eqtl_pred_model_df <- table_filt_df %>%
    filter(sig_eqtl_any == T) %>%
    pull(gene_id.v26) %>%
    get_prediction_model_rsids(., tiss_filter=tiss_filter, "eqtl") %>%
    left_join(meta_ukb_df, by=c("rsid" = "rsID")) %>%
    left_join(pooled_ukb_df, by=c("rsid" = "rsID")) %>%
    mutate(beta_diff = abs(beta_meta - beta_pool),
           SE_diff = abs(SE_meta - SE_pool),
           p_diff = abs(p_meta - p_pool)) %>%
    filter(REF_meta == REF_pool & ALT_meta == ALT_pool)

  sqtl_pred_model_df <- table_filt_df %>%
    filter(sig_sqtl_any == T) %>%
    pull(gene_id.v26) %>%
    get_prediction_model_rsids(., tiss_filter=tiss_filter, "sqtl") %>%
    left_join(meta_ukb_df, by=c("rsid" = "rsID")) %>%
    left_join(pooled_ukb_df, by=c("rsid" = "rsID")) %>%
    mutate(beta_diff = abs(beta_meta - beta_pool),
           SE_diff = abs(SE_meta - SE_pool),
           p_diff = abs(p_meta - p_pool)) %>%
    filter(REF_meta == REF_pool & ALT_meta == ALT_pool)

  sqtl_model_num_uniq_vars <- sqtl_pred_model_df %>% pull(varID) %>% unique() %>% length()
  eqtl_model_num_uniq_vars <- eqtl_pred_model_df %>% pull(varID) %>% unique() %>% length()

  sqtl_model_num_uniq_matching_vars <- sqtl_pred_model_df %>%
    dplyr::filter(!is.na(beta_diff)) %>%
    pull(varID) %>%
    unique() %>% 
    length()

  eqtl_model_num_uniq_matching_vars <- eqtl_pred_model_df %>%
    dplyr::filter(!is.na(beta_diff)) %>%
    pull(varID) %>%
    unique() %>% 
    length()

  eqtl_pred_model_matching_df <- eqtl_pred_model_df %>%
    dplyr::filter(!is.na(beta_diff)) %>%
    rename(gene_id.v26 = gene) %>%
    left_join(names_filt_df, by="gene_id.v26") %>%
    select(all_of(c("gene_id.v26", "gene_name.v26", "gene_name.v40", "varID",
    "rsid", "by_hand_rsid", "tissue", "weight", "beta_meta", "SE_meta", "p_meta",
    "beta_pool","SE_pool", "p_pool", "beta_diff", "SE_diff", "p_diff")))

  sqtl_pred_model_matching_df <- sqtl_pred_model_df %>%
    dplyr::filter(!rlang::are_na(beta_diff)) %>%
    rename(gene_id.v26 = gene_id) %>%
    left_join(names_filt_df, by="gene_id.v26") %>%
    select(all_of(c("gene_id.v26", "gene_name.v26", "gene_name.v40", "varID",
    "rsid", "by_hand_rsid", "tissue", "weight", "beta_meta", "SE_meta", "p_meta",
    "beta_pool","SE_pool", "p_pool", "beta_diff", "SE_diff", "p_diff")))
  
  eqtl_uniq_df <- eqtl_pred_model_matching_df %>%
    select(all_of(c("varID", "rsid", "by_hand_rsid", "beta_meta", "SE_meta", "p_meta",
    "beta_pool","SE_pool", "p_pool",
     "beta_diff", "SE_diff", "p_diff"))) %>%
    distinct() %>%
    arrange(varID)

  sqtl_uniq_df <- sqtl_pred_model_matching_df %>%
    select(all_of(c("varID", "rsid", "by_hand_rsid", "beta_meta", "SE_meta", "p_meta",
    "beta_pool","SE_pool", "p_pool",
     "beta_diff", "SE_diff", "p_diff"))) %>%
    distinct() %>%
    arrange(varID)

  combined_pred_model_matching_df <- bind_rows(eqtl_pred_model_matching_df, sqtl_pred_model_matching_df) %>%
    distinct() %>%
    arrange(varID)

  combined_uniq_df <- bind_rows(eqtl_uniq_df, sqtl_uniq_df) %>%
    distinct() %>%
    arrange(varID)

  write_delim(eqtl_pred_model_matching_df, glue("{odir}all_eqtl_pred_model_matching.tsv"), delim="\t")
  write_delim(sqtl_pred_model_matching_df, glue("{odir}all_sqtl_pred_model_matching.tsv"), delim="\t")
  write_delim(combined_pred_model_matching_df, glue("{odir}all_pred_model_matching.tsv"), delim="\t")

  write_delim(eqtl_uniq_df, glue("{odir}uniq_var_eqtl_model_matching.tsv"), delim="\t")
  write_delim(sqtl_uniq_df, glue("{odir}uniq_var_sqtl_model_matching.tsv"), delim="\t")
  write_delim(combined_uniq_df, glue("{odir}uniq_var_pred_model_matching.tsv"), delim="\t")
}

get_prediction_model_rsids <- function(gene_ids.v26, tiss_filter, eqtl_or_sqtl){
  require(tidyverse)
  require(rlang)
  print(glue("Getting {eqtl_or_sqtl} prediction model rsids for {length(gene_ids.v26)} genes"))
  want_tissues <- fread(tiss_filter, header=F) %>% pull(V1)
  rsid_fix_df <- fread("../input/hg38_varid_to_rsid.csv") # Manually curated file to map non-rsid entries to rsids

  if (eqtl_or_sqtl == "sqtl"){
    pheno_df <- fread("../input/combine_phenotype_groups.txt.gz") %>%
      filter(gene_id %in% gene_ids.v26) %>%
      separate(intron_id, c("tissue", "intron_id"), sep="\\.", remove=T)
  }

  extract_from_tiss_db <- function(want_tiss){
    require(tidyverse)
    require(rlang)
    print(want_tiss)
    con <- DBI::dbConnect(RSQLite::SQLite(), glue("../input/gtex_v8_{eqtl_or_sqtl}_dbs_mashr/mashr_{want_tiss}.db"))

    weights_db <- tbl(con, "weights")

    if (eqtl_or_sqtl == "eqtl"){
      weights_filt_df <- weights_db %>%
        filter(gene %in% gene_ids.v26) %>%
        as_tibble() %>%
        mutate(tissue = want_tiss)
    } else if (eqtl_or_sqtl == "sqtl"){
      pheno_tiss_df <- pheno_df %>%
        filter(tissue == want_tiss)
      weights_filt_df <- weights_db %>%
        as_tibble() %>%
        rename(intron_id = gene) %>%
        inner_join(pheno_tiss_df, by="intron_id") %>%
        filter(gene_id %in% gene_ids.v26)
    }
    DBI::dbDisconnect(con)
    
    # The publicly available gtex v8 weights databases have some non-rsid
    # entries in the `rsid` column, but looking at
    # https://www.ncbi.nlm.nih.gov/snp/, the rsids can actually be found.
    #
    # Thus, I manually curate a file to map these non-rsid entries to rsids
    non_rsid_df <- weights_filt_df %>% 
      filter(!grepl("rs[0-9]+", rsid)) %>%
      rename(not_rsid = rsid) %>%
      left_join(rsid_fix_df, by="not_rsid") %>% # Use the rsid_fix_df to map non-rsid entries to rsids
      select(-all_of(c("not_rsid")))

    rsid_df <- weights_filt_df %>% 
      filter(grepl("rs[0-9]+", rsid))

    still_no_rsid <- non_rsid_df %>%
      filter(rlang::are_na(rsid))

    if (nrow(still_no_rsid) > 0){
      print(still_no_rsid$varID %>% unique())
      stop("There are still unmapped non-rsid entries in the weights database")
    } else {
      non_rsid_df <- non_rsid_df %>%
        mutate(by_hand_rsid = TRUE)
      rsid_df <- rsid_df %>%
        mutate(by_hand_rsid = FALSE)
      
      return(rbind(non_rsid_df, rsid_df))
    }
  }
  return_df <- want_tissues %>% map_dfr(~extract_from_tiss_db(.))
  return(return_df)
}



if (!interactive()) {
  main()
} else {
  print("Sourced interactively.")
  main()
}