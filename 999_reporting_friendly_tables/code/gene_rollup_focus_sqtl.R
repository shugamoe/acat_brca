FOCUS_RESULTS <- list.files(path="../../05_focus/output/", pattern="focus_sqtl_11tiss_intrinsic_subtype_.*.tsv",
                            full.names=TRUE)
# FOCUS_RESULTS <- list.files(path="../../05_focus/output/", pattern="focus_sqtl_11tiss_meta_analysis_BCAC_UKB_ovr_brca_no_pval_filt_by_tissue.focus.tsv",
#                             full.names=TRUE)

library(tidyverse)
library(data.table)
library(glue)
GROUPINGS <- fread("../input/combine_phenotype_groups.txt.gz") %>% 
  as_tibble() %>% 
  rename(group = gene_id)
GROUPINGS[c("tissue", "intron")] <- stringr::str_split_fixed(GROUPINGS$intron_id, pattern="\\.", 2)
GROUPINGS$tissue <- tolower(GROUPINGS$tissue)

# Guimin followups
# GM_FUP <- c("../../05_focus/output//focus_sqtl_11tiss_meta_analysis_BCAC_UKB_ovr_brca_by_tissue.focus.tsv", 
# "../../05_focus/output//focus_sqtl_11tiss_meta_analysis_BCAC_UKB_ovr_brca.focus.tsv")
GM_FUP <- c("../../05_focus/output//focus_sqtl_11tiss_meta_analysis_BCAC_UKB_ovr_brca.focus.tsv")
GM_FUP_GENES <- c(
# NULL.MODEL pip?
"ENSG00000146833.15",
"ENSG00000176402.5",
"ENSG00000160862.12",
"ENSG00000078319.9",
"ENSG00000242294.6",
"ENSG00000121716.20",
"ENSG00000085514.15",
"ENSG00000078487.17",
"ENSG00000166925.8",
"ENSG00000166924.8",

# in same LD?
"ENSG00000164440.14",
"ENSG00000226571.1",

# same LD?
"ENSG00000116288.12",
"ENSG00000116285.12")

main <- function(input_focus_sqtl_path,
                 odir="../output",
                 breast_filt=F){
  require(tidyverse)
  require(data.table)
  require(glue)
  require(purrr)

  print(input_focus_sqtl_path)
  fname_only <- basename(input_focus_sqtl_path)

  if (input_focus_sqtl_path %in% GM_FUP){
    followup_file <- TRUE
  } else {
    followup_file <- FALSE
  }

  # ens_gene_id ens_tx_id   mol_name    tissue  ref_name    type    chrom   tx_start    tx_stop inference   cv.R2   cv.R2.pval  twas_z  pip in_cred_set region
  focus_sqtl_df <- fread(input_focus_sqtl_path)

  get_output_cols <- function(tissue_target=NA){
    if (!is.na(tissue_target)){
      focus_sqtl_df_filt <- focus_sqtl_df %>% filter(tissue_target == !!tissue_target)
    } else {
      focus_sqtl_df_filt <- focus_sqtl_df
    }

    if (followup_file == TRUE){
      focus_fup_grouped <- inner_join(GROUPINGS %>%
      filter(group %in% GM_FUP_GENES), focus_sqtl_df_filt %>% rename(intron=ens_gene_id), by = c("intron", "tissue")) %>%
        select(-c("intron_id", "gtex_intron_id")) %>%
        select(any_of(c("group", "intron", "tissue", "region", "tissue_target", "pip")))
      focus_sqtl_df_null <- focus_sqtl_df_filt %>%
        rename(intron = ens_gene_id,
               null_model_tissue_region_pip = pip) %>%
        select(-any_of(c("tissue", "tissue_target"))) %>% # tissue is '' for NULL.MODEL and tissue_target is irrelevant
        filter(intron == "NULL.MODEL",
        region %in% focus_fup_grouped$region)

      rdf_guimin <- inner_join(focus_fup_grouped, focus_sqtl_df_null, by = c("region")) %>%
        select(any_of(c("group", "intron", "tissue", "region", "pip", "null_model_tissue_region_pip", "tissue_target"))) %>%
        arrange(region, group, tissue)
      
      oname <- glue("../output/guimin_tables/focus_sqtl_followup_pip_ld_blocks_tissue_target_{tissue_target}_breast_filt_{breast_filt}.tsv")
      write_delim(rdf_guimin, oname, delim="\t")
    }


    focus_sqtl_df_filt <- focus_sqtl_df_filt %>%
      filter(ens_gene_id != "NULL.MODEL") %>% # Leave out NULL.MODEL, investigate raw FOCUS results if curious
      select(all_of(c("ens_gene_id", "chrom", "twas_z", "pip", "in_cred_set", "tissue"))) %>%
      rename(intron = ens_gene_id)
  
    if (isTRUE(breast_filt)){
      focus_sqtl_df_filt <- focus_sqtl_df_filt %>% filter(tissue == "breast_mammary_tissue")
    }

    focus_grouped <- inner_join(GROUPINGS, focus_sqtl_df_filt, by = c("intron", "tissue")) %>%
      select(-c("intron_id", "gtex_intron_id"))

    focus_highest_intron_pip <- focus_grouped %>% 
      group_by(group) %>%
      summarize(max_pip = max(pip, na.rm=T))

    focus_highest_intron_pip_tissue <- focus_grouped %>%
      left_join(focus_highest_intron_pip, by = c("group")) %>%
      filter(pip == max_pip) %>%
      rename(max_pip_intron = intron) %>%
      group_by(group, max_pip) %>%
      summarize(max_pip_introns = paste(unique(max_pip_intron) %>% sort(), collapse = ", "),
                max_pip_tissues = paste(unique(tissue) %>% sort(), collapse = ", "),
                max_pip_num_in_cred_set = sum(in_cred_set, na.rm=T))

    print(names(focus_highest_intron_pip_tissue))

    # Remove duplicates from FOCUS (genes spanning multiple LD regions)
    # print(names(focus_sqtl_df_filt))
    dup_key <- focus_grouped %>%
      group_by(group) %>%
      summarize(introns_in_focus = n(),
                introns_in_cred_set = sum(in_cred_set, na.rm=T),
                highest_pip = max(pip, na.rm=T),
                highest_abs_twas_z = max(abs(twas_z), na.rm=T)) %>%
      mutate(introns_in_focus = ifelse(highest_pip == -Inf, 0, introns_in_focus),
            highest_pip = ifelse(highest_pip == -Inf, NA, highest_pip),
            introns_in_cred_set = ifelse(highest_pip == -Inf, 0, introns_in_cred_set),
            highest_abs_twas_z = ifelse(highest_abs_twas_z == -Inf, NA, highest_abs_twas_z),
            
            )

    kill_cols <- c("twas_z", "pip", "intron", "tissue", "cluster_id", "in_cred_set", "max_pip")
    t_out_df <- left_join(focus_grouped, dup_key, by=c("group")) %>%
      left_join(focus_highest_intron_pip_tissue) %>%
      rename(ens_gene_id = group) %>%
      select(-any_of(kill_cols)) %>%
      distinct()
    
    # if (is.na(tissue_target)){
    #   browser()
    # }
    # browser()
    if (!is.na(tissue_target)){
      t_out_df <- t_out_df %>% 
        rename(!!glue("{tissue_target}_introns_in_focus") := introns_in_focus,
               !!glue("{tissue_target}_introns_in_cred_set") := introns_in_cred_set,
               !!glue("{tissue_target}_highest_pip") := highest_pip,
               !!glue("{tissue_target}_highest_abs_twas_z") := highest_abs_twas_z,
               !!glue("{tissue_target}_max_pip_introns") := max_pip_introns,
               !!glue("{tissue_target}_max_pip_tissues") := max_pip_tissues
        )
      return(t_out_df %>% select(-all_of(c(glue("{tissue_target}_max_pip_tissues"), "max_pip_num_in_cred_set"))))
    } else {
      browser()
      return(t_out_df)
    }
  }

  if (isTRUE(breast_filt)){
    breast_tag <- "breast_only_"
  } else {
    breast_tag <- ""
  }

  # Check what type of focus file, by tissue or individ/ual run?
  if ("tissue_target" %in% names(focus_sqtl_df)){
    if (breast_filt == T){
      return(NULL) # Don't need to make breast tissue only files for by tissue files
    }

    tissues_present <- focus_sqtl_df %>%
      pull(tissue_target) %>%
      unique()
    tissues_present <- c(NA, tissues_present) # Add NA to get max pip tissues across all tissue-exclusive models

    df_list <- tissues_present %>%
      map(~get_output_cols(tissue_target=.))

    t_out_df <- df_list %>% reduce(full_join, by = c("ens_gene_id", "chrom"))
  } else {
    t_out_df <- get_output_cols()
  }
  # browser()
 
  opath <- file.path(odir, glue("gene_rolled_{breast_tag}{fname_only}"))
  write_delim(t_out_df, opath, delim="\t")
}

if (!interactive()){
  library(purrr)
  FOCUS_RESULTS %>% walk(~main(.))
  # FOCUS_RESULTS %>% walk(~main(., breast_filt=T))
} else {
  library(purrr)
  FOCUS_RESULTS %>% walk(~main(.))
  # gimme <- main(FOCUS_RESULTS[7])
  # wantthis <- main(FOCUS_RESULTS[8])
}
