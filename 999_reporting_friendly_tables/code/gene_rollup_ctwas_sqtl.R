# FOCUS_RESULTS <- list.files(path="../../06_focus/output/", pattern="focus_sqtl.*.focus.tsv",
#                             full.names=TRUE)

library(tidyverse)
library(data.table)
library(glue)
GROUPINGS <- fread("../input/combine_phenotype_groups.txt.gz") %>% 
  as_tibble() %>% 
  rename(group = gene_id)
GROUPINGS[c("tissue", "intron")] <- stringr::str_split_fixed(GROUPINGS$intron_id, pattern="\\.", 2)

main <- function(study_short,
                 tiss_list="../input/tissue_lists/11tiss.txt",
                 odir="../output",
                 breast_filt=F,
                 by_chrom=F,
                 test=T){
  require(tidyverse)
  require(data.table)
  require(glue)
  require(purrr)

  study_long <- study_transform(study_short, to="long")

  want_tissues <- fread(tiss_list, header=F)[["V1"]]
  if (test){
    want_tissues <- want_tissues[c(1,2)]
  }

  ctwas_intron_df <- tibble()
  print(study_long)
  for (cur_tiss in want_tissues){
    if (test){
      fpath <- "/gpfs/data/biostat-lab/Gao_group/Julian/AD_PWAS_custom_models/02_ctwas/output/ctwas_rss/schwartz_gwas_wingo_nc_2022_fusion_best.susieIrss.txt"
      offset <- runif(1)
    } else {
      fpath <- glue("../../07_ctwas/output/ctwas_rss/{study_short}__{cur_tiss}__sqtl.susieIrss.txt")
      if (by_chrom){
      } else {
        fpath <- glue("../../07_ctwas/output/ctwas_rss/{study_short}__{cur_tiss}__sqtl.susieIrss.txt")
        if (!file.exists(fpath)){
          print(glue("missing {fpath}"))
          next
        }
        offset <- 0
      }
    }
    if (by_chrom){
      all_chroms_df <- tibble()
      want_files <- list.files(path="../../07_ctwas/output/ctwas_rss/", pattern=glue("{study_short}__{cur_tiss}__sqtl_chr\\d\\d*.susieIrss.txt"), full.names=T)
      for (cur_file in want_files){
        cur_df <- fread(cur_file) %>%
          filter(type == "gene") %>%
          mutate(tissue = cur_tiss) %>%
          select(all_of(c("chrom", "id", "tissue", "susie_pip", "region_tag2")))
        all_chroms_df <- bind_rows(all_chroms_df, cur_df)
      }
      ctwas_intron_df <- bind_rows(ctwas_intron_df, all_chroms_df)
    } else {
      cur_df <- fread(fpath) %>%
          filter(type == "gene")
      cur_df <- cur_df %>%
        mutate(tissue = cur_tiss,
               susie_pip = susie_pip + offset) %>%
        select(all_of(c("chrom", "id", "tissue", "susie_pip", "region_tag2")))
      ctwas_intron_df <- bind_rows(ctwas_intron_df, cur_df)
    }
  }
  if (test){
    t_introns <- c("intron_10_100003023_100009838",
    "intron_10_1001013_1005818",
    "intron_10_100164095_100174208",
    "intron_10_100190968_100193298",
    "intron_10_100192693_100193298",
    "intron_10_100245913_100246795",
    "intron_10_100246130_100246795",
    "intron_10_100246935_100253421",
    "intron_10_100253539_100256262",
    "intron_10_100253539_100260218")
    nr <- nrow(ctwas_intron_df)
    set.seed(117)
    ctwas_intron_df[["id"]] <- sample(t_introns, nr, replace=T)
    ctwas_intron_df <- ctwas_intron_df[c(1:6),]
  }
  fname_only <- basename(fpath)

  get_output_cols <- function(tissue_target=NA){
    if (!is.na(tissue_target)){
      ctwas_intron_df_filt <- ctwas_intron_df %>% filter(tissue == !!tissue_target)
    } else {
      ctwas_intron_df_filt <- ctwas_intron_df
    }

    ctwas_intron_df_filt <- ctwas_intron_df_filt %>%
      select(all_of(c("id", "chrom", "susie_pip", "tissue", "region_tag2"))) %>%
      rename(intron = id)
  
    if (isTRUE(breast_filt)){
      ctwas_intron_df_filt <- ctwas_intron_df_filt %>% filter(tissue == "Breast_Mammary_Tissue")
    }

    ctwas_grouped <- inner_join(GROUPINGS, ctwas_intron_df_filt, by = c("intron", "tissue")) %>%
      select(-c("intron_id", "gtex_intron_id"))

    ctwas_highest_intron_susie_pip <- ctwas_grouped %>% 
      group_by(group) %>%
      summarize(max_susie_pip = max(susie_pip, na.rm=T),
                max_region_tag = max(region_tag2, na.rm=T),
                min_region_tag = min(region_tag2, na.rm=T),
                )

    ctwas_highest_intron_susie_pip_tissue <- ctwas_grouped %>%
      left_join(ctwas_highest_intron_susie_pip, by = c("group")) %>%
      filter(susie_pip == max_susie_pip) %>%
      rename(max_susie_pip_intron = intron) %>%
      group_by(group, max_susie_pip) %>%
      summarize(max_susie_pip_introns = paste(unique(max_susie_pip_intron) %>% sort(), collapse = ", "),
                highest_susie_pip_tissues = paste(unique(tissue) %>% sort(), collapse = ", "))

    print(names(ctwas_highest_intron_susie_pip_tissue))

    # Remove duplicates from FOCUS (genes spanning multiple LD regions)
    # print(names(ctwas_intron_df_filt))
    dup_key <- ctwas_grouped %>%
      group_by(group) %>%
      summarize(highest_susie_pip = max(susie_pip, na.rm=T)) %>%
      mutate(highest_susie_pip = ifelse(highest_susie_pip == -Inf, NA, highest_susie_pip))

    kill_cols <- c("susie_pip", "intron", "tissue", "max_susie_pip", "cluster_id", "region_tag2")
    t_out_df <- left_join(ctwas_grouped, dup_key, by=c("group")) %>%
      left_join(ctwas_highest_intron_susie_pip_tissue) %>%
      rename(ens_gene_id = group) %>%
      select(-any_of(kill_cols)) %>%
      distinct()
    
    if (!is.na(tissue_target)){
      t_out_df <- t_out_df %>% 
        rename(
               !!glue("{tissue_target}_highest_susie_pip") := highest_susie_pip,
               !!glue("{tissue_target}_max_susie_pip_introns") := max_susie_pip_introns,
               !!glue("{tissue_target}_highest_susie_pip_tissues") := highest_susie_pip_tissues,
               !!glue("{tissue_target}_max_region_tag") := max_region_tag,
               !!glue("{tissue_target}_min_region_tag") := min_region_tag
        )
      return(t_out_df %>% select(-all_of(c(glue("{tissue_target}_highest_susie_pip_tissues")))))
    } else {
      return(t_out_df)
    }
  }

  if (isTRUE(breast_filt)){
    breast_tag <- "breast_only_"
  } else {
    breast_tag <- ""
  }

  if (breast_filt == T){
    return(NULL) # Don't need to make breast tissue only files for by tissue files
  }

  tissues_present <- ctwas_intron_df %>%
    pull(tissue) %>%
    unique()
  tissues_present <- c(NA, tissues_present) # Add NA to get max pip tissues across all tissue-exclusive models

  df_list <- tissues_present %>%
    map(~get_output_cols(tissue_target=.))

  t_out_df <- df_list %>% reduce(full_join, by = c("ens_gene_id", "chrom"))
 
  fname_only <- fname_only %>%
    str_replace(., cur_tiss, "") %>%
    str_replace(., study_short, study_long) %>%
    str_replace(., "____", "_")
  opath <- file.path(odir, glue("gene_rolled_ctwas_{breast_tag}{fname_only}"))
  # browser()
  write_delim(t_out_df, opath, delim="\t")
}

study_transform <- function(study_name, to){
  short_to_long <- switch(study_name, 
                      "bcac_erp" = "BCAC_ERPOS_BreastCancer_EUR",
                      "meta_bcac_cimba_erneg" = "meta_analysis_BCAC_CIMBA_erneg_brca",
                      "bcac_ovr" = "BCAC_Overall_BreastCancer_EUR",
                      "meta_bcac_ukb_ovr" = "meta_analysis_BCAC_UKB_ovr_brca",
                      # "short_to_long_not_found"
                      study_name
                        )
  long_to_short <- switch(study_name, 
                      "BCAC_ERPOS_BreastCancer_EUR" = "bcac_erp",
                      "meta_analysis_BCAC_CIMBA_erneg_brca" = "meta_bcac_cimba_erneg",
                      "BCAC_Overall_BreastCancer_EUR" = "bcac_ovr",
                      "meta_analysis_BCAC_UKB_ovr_brca" = "meta_bcac_ukb_ovr",
                      "long_to_short_not_found"
                        )
  if (to == "short"){
    return(long_to_short)
  } else if (to == "long"){
    return(short_to_long)
  }
}

get_intron_tiss_from_ensg <- function(ensg_v26_ids,
                                      tiss_list="../input/tissue_lists/11tiss.txt",
                                      grouping=GROUPINGS){
  want_tissues <- fread(tiss_list, header=F)[["V1"]]

  want_grouping <- grouping %>%
    filter(group %in% ensg_v26_ids,
           tissue %in% want_tissues)
  browser()

}

if (!interactive()){
  library(purrr)
  paste0("intrinsic_subtype_", 1:5) %>% walk(~main(., test=F, by_chrom=T))
  main(study_short="meta_bcac_ukb_ovr", test=F, by_chrom=T)
} else {
  main(study_short="meta_bcac_ukb_ovr", test=F, by_chrom=T)
  # get_intron_tiss_from_ensg(c("ENSG00000197935.6 ", "ENSG00000226067.6"))
}
