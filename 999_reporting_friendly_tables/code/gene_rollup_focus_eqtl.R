# FOCUS_RESULTS <- list.files(path="../../05_focus/output/", pattern="focus_eqtl_11tiss_.*.focus.tsv",
#                             full.names=TRUE)
FOCUS_RESULTS <- list.files(path="../../05_focus/output/", pattern="focus_eqtl_11tiss_meta_analysis_BCAC_UKB_ovr_brca_no_pval_filt_by_tissue.focus.tsv",
                            full.names=TRUE)

main <- function(input_focus_eqtl_path,
                 odir="../output",
                 breast_filt=F){
  require(tidyverse)
  require(data.table)
  require(glue)

  print(input_focus_eqtl_path)
  fname_only <- basename(input_focus_eqtl_path)

  focus_eqtl_df <- fread(input_focus_eqtl_path)
  
  get_output_cols <- function(tissue_target=NA){
    if (!is.na(tissue_target)){
      focus_eqtl_df_filt <- focus_eqtl_df %>% filter(tissue_target == !!tissue_target)
    } else {
      focus_eqtl_df_filt <- focus_eqtl_df
    }

    focus_eqtl_df_filt <- focus_eqtl_df_filt %>%
      filter(ens_gene_id != "NULL.MODEL") %>% # Leave out NULL.MODEL, investigate raw FOCUS results if curious %>%
      mutate(ens_gene_id_no_decimal = str_extract(ens_gene_id, "(ENSG\\d{1,})"))

    if (isTRUE(breast_filt)){
      focus_eqtl_df_filt <- focus_eqtl_df_filt %>% filter(tissue == "breast_mammary_tissue")
    }

    # Remove duplicates from FOCUS (genes spanning multiple LD regions)
    # print(names(focus_eqtl_df_filt))
    dup_key <- focus_eqtl_df_filt %>%
      group_by(ens_gene_id_no_decimal, chrom) %>%
      summarize(num_in_focus = n(),
                num_in_cred_set = sum(in_cred_set, na.rm=T),
                highest_pip = max(pip, na.rm=T),
                highest_abs_twas_z = max(abs(twas_z), na.rm=T)) %>%
      mutate(num_in_focus = ifelse(highest_pip == -Inf, 0, num_in_focus),
            highest_pip = ifelse(highest_pip == -Inf, NA, highest_pip),
            num_in_cred_set = ifelse(highest_pip == -Inf, 0, num_in_cred_set),
            highest_abs_twas_z = ifelse(highest_abs_twas_z == -Inf, NA, highest_abs_twas_z))

    # These columns have been grouped above, or leaving them in would cause duplicates
    kill_cols <- c("twas_z", "pip", "in_cred_set", "region", "tissue",
    "tx_start", "tx_stop", "cv.R2", "cv.R2.pval", "ens_tx_id", "tissue_target", "inference", "ref_name")

    t_out_df <- focus_eqtl_df_filt %>%
      left_join(dup_key, by = c("ens_gene_id_no_decimal", "chrom"))

    if (is.na(tissue_target) & isFALSE(breast_filt)){
      eqtl_max_pip <- focus_eqtl_df_filt %>%
        group_by(ens_gene_id_no_decimal) %>%
        summarize(max_pip = max(pip, na.rm=T))
      tiss_max_pip = focus_eqtl_df_filt %>%
        left_join(eqtl_max_pip, by = c("ens_gene_id_no_decimal")) %>%
        filter(pip == max_pip) %>%
        group_by(ens_gene_id_no_decimal, max_pip) %>%
        summarize(highest_pip_tissues = paste(unique(tissue) %>% sort(), collapse = ", ")) %>%
        select(-all_of(c("max_pip")))

      t_out_df <- t_out_df %>%
        left_join(tiss_max_pip, by = c("ens_gene_id_no_decimal"))
    }

    t_out_df <- t_out_df %>%
      select(-any_of(kill_cols)) %>%
      distinct()

    if (!is.na(tissue_target)){
      t_out_df <- t_out_df %>% 
        rename(!!glue("{tissue_target}_num_in_focus") := num_in_focus,
               !!glue("{tissue_target}_num_in_cred_set") := num_in_cred_set,
               !!glue("{tissue_target}_highest_pip") := highest_pip,
               !!glue("{tissue_target}_highest_abs_twas_z") := highest_abs_twas_z
        )
      return(t_out_df)
    } else {
      return(t_out_df)
    }
  }

  if (isTRUE(breast_filt)){
    breast_tag <- "breast_only_"
  } else {
    breast_tag <- ""
  }

  # Check what type of focus file, by tissue or individual run?
  if ("tissue_target" %in% names(focus_eqtl_df)){
    if (breast_filt == T){
      return(NULL) # Don't need to make breast tissue only files for by tissue files
    }

    tissues_present <- focus_eqtl_df %>%
      pull(tissue_target) %>%
      unique()
    tissues_present <- c(NA, tissues_present) # Add NA to get max pip tissues across all tissue-exclusive models

    df_list <- tissues_present %>%
      map(~get_output_cols(tissue_target=.))

    t_out_df <- df_list %>% reduce(full_join, by = c("ens_gene_id_no_decimal",
    "chrom", "ens_gene_id", "mol_name", "type"))
  } else {
    t_out_df <- get_output_cols()
  }
  #browser()
  opath <- file.path(odir, glue("gene_rolled_{breast_tag}{fname_only}"))
  write_delim(t_out_df, opath, delim="\t")
  # return(t_out_df)
}

if (!interactive()) {
  library(purrr)
  print(FOCUS_RESULTS)
  FOCUS_RESULTS %>% walk(~main(.))
  # FOCUS_RESULTS %>% walk(~main(., breast_filt=T))
} else {
  # gimme <- main(FOCUS_RESULTS[1])
  library(purrr)
  FOCUS_RESULTS %>% walk(~main(.))
}
