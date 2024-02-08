main <- function(study_short,
                 # tiss_list="../input/tissue_lists/11tiss.txt",
                 tiss_list="../input/tissue_lists/Breast_Mammary.txt",
                 odir="../output",
                 breast_filt=F,
                 by_chrom=F,
								 z_cutoff=NULL,
                 test=T){
  require(tidyverse)
  require(data.table)
  require(glue)

  study_long <- study_transform(study_short, to="long")

  want_tissues <- fread(tiss_list, header=F)[["V1"]]
  if (test){
    want_tissues <- want_tissues[c(1,2)]
  }

  ctwas_gene_df <- tibble()
  print(study_long)
  for (cur_tiss in want_tissues){
    if (test){
      fpath <- "/gpfs/data/biostat-lab/Gao_group/Julian/AD_PWAS_custom_models/02_ctwas/output/ctwas_rss/schwartz_gwas_wingo_nc_2022_fusion_best.susieIrss.txt"
      offset <- rnorm(1)
    } else {
      fpath <- glue("../../07_ctwas/output/ctwas_rss/{study_short}__max_mag_snp_z{z_cutoff}__{cur_tiss}__eqtl.susieIrss.txt")
      if (by_chrom){
        print(glue("Getting study by chromosome tissue: {cur_tiss}"))
      } else {
        fpath <- glue("../../07_ctwas/output/ctwas_rss/{study_short}__max_mag_snp_z{z_cutoff}__{cur_tiss}__eqtl.susieIrss.txt")
        if (!file.exists(fpath)){
          print(glue("missing {fpath}"))
          next
        }
        offset <- 0
      }
    }
    if (by_chrom){
      all_chroms_df <- tibble()
      want_files <- list.files(path="../../07_ctwas/output/ctwas_rss/", pattern=glue("{study_short}__max_mag_snp_z{z_cutoff}__{cur_tiss}__eqtl_chr\\d\\d*.susieIrss.txt"),
                               full.names=T)
      for (cur_file in want_files){
        cur_df <- fread(cur_file) %>%
          filter(type == "gene") %>%
          mutate(tissue_target = cur_tiss)
        all_chroms_df <- bind_rows(all_chroms_df, cur_df)
      }
      ctwas_gene_df <- bind_rows(ctwas_gene_df, all_chroms_df)
    } else {
      cur_df <- fread(fpath) %>%
          filter(type == "gene")
      cur_df <- cur_df %>%
        mutate(tissue_target = cur_tiss,
               susie_pip = susie_pip + offset) %>%
        select(all_of(c("chrom", "id", "tissue_target", "susie_pip", "region_tag2")))
      ctwas_gene_df <- bind_rows(ctwas_gene_df, cur_df)
    }
  }
  fname_only <- basename(fpath)

  get_output_cols <- function(tissue_target=NA){
    print(glue("Get output columns by {tissue_target}"))
    if (!is.na(tissue_target)){
      ctwas_gene_df_filt <- ctwas_gene_df %>% filter(tissue_target == !!tissue_target)
    } else {
      ctwas_gene_df_filt <- ctwas_gene_df
    }

    if (isTRUE(breast_filt)){
      ctwas_gene_df_filt <- ctwas_gene_df_filt %>% filter(tissue_target == "Breast_Mammary_Tissue")
    }

    # Remove duplicates from FOCUS (genes spanning multiple LD regions)
    # print(names(ctwas_gene_df_filt))
    dup_key <- ctwas_gene_df_filt %>%
      group_by(id, chrom) %>%
      summarize(highest_susie_pip = max(susie_pip, na.rm=T),
               min_region_tag = min(region_tag2, na.rm=T),
               max_region_tag = max(region_tag2, na.rm=T)
                ) %>%
      mutate(highest_susie_pip = ifelse(highest_susie_pip == -Inf, NA, highest_susie_pip))

    # These columns have been grouped above, or leaving them in would cause duplicates
    kill_cols <- c("susie_pip", "type", "pos", "tissue_target", "region_tag1", "region_tag2", "cs_index", "mu2")

    t_out_df <- ctwas_gene_df_filt %>%
      left_join(dup_key, by = c("id", "chrom"))

    if (is.na(tissue_target) & isFALSE(breast_filt)){
      eqtl_max_pip <- ctwas_gene_df_filt %>%
        group_by(id) %>%
        summarize(max_susie_pip = max(susie_pip, na.rm=T))
      tiss_max_pip = ctwas_gene_df_filt %>%
        left_join(eqtl_max_pip, by = c("id")) %>%
        filter(susie_pip == max_susie_pip) %>%
        group_by(id, max_susie_pip) %>%
        summarize(highest_susie_pip_tissues = paste(unique(tissue_target) %>% sort(), collapse = ", ")) %>%
        select(-all_of(c("max_susie_pip")))

      t_out_df <- t_out_df %>%
        left_join(tiss_max_pip, by = c("id"))
    }

    t_out_df <- t_out_df %>%
      select(-any_of(kill_cols)) %>%
      distinct()

    if (!is.na(tissue_target)){
      t_out_df <- t_out_df %>% 
        rename(!!glue("{tissue_target}_highest_susie_pip") := highest_susie_pip,
               !!glue("{tissue_target}_max_region_tag") := max_region_tag,
               !!glue("{tissue_target}_min_region_tag") := min_region_tag
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
  if (breast_filt == T){
    return(NULL) # Don't need to make breast tissue only files for by tissue files
  }

  tissues_present <- ctwas_gene_df %>%
    pull(tissue_target) %>%
    unique()
  tissues_present <- c(NA, tissues_present) # Add NA to get max pip tissues across all tissue-exclusive models

  df_list <- tissues_present %>%
    map(~get_output_cols(tissue_target=.))

  t_out_df <- df_list %>% reduce(full_join, by = c("id",
  "chrom"))

  fname_only <- fname_only %>%
    str_replace(., cur_tiss, "") %>%
    str_replace(., study_short, study_long) %>%
    str_replace(., "____", "_")
   
  opath <- file.path(odir, glue("gene_rolled_ctwas_{breast_tag}{fname_only}"))
  # browser()
  write_delim(t_out_df, opath, delim="\t")
  # return(t_out_df)
}

study_transform <- function(study_name, to){
  short_to_long <- switch(study_name, 
                      "bcac_erp" = "BCAC_ERPOS_BreastCancer_EUR",
                      "meta_bcac_cimba_erneg" = "meta_analysis_BCAC_CIMBA_erneg_brca",
                      "bcac_ovr" = "BCAC_Overall_BreastCancer_EUR",
                      "meta_bcac_ukb_ovr" = "meta_analysis_BCAC_UKB_ovr_brca",
                      "meta_bcac_ukb_ovr_no_pval_filt" = "meta_analysis_BCAC_UKB_ovr_brca_no_pval_filt",
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

TEST_CUTOFFS <- c(5, 10, 15, 20, 25, 30, 500)
if (!interactive()) {
  library(purrr)
  # FOCUS_RESULTS %>% walk(~main(.))
  # FOCUS_RESULTS %>% walk(~main(., breast_filt=T))
  # paste0("intrinsic_subtype_", 1:5) %>% walk(~main(., test=F, by_chrom=T))
  # main(study_short="meta_bcac_ukb_ovr", test=F, by_chrom=T)
  # main(study_short="intrinsic_subtype_1", test=F, by_chrom=T)
} else {
  library(purrr)
	for (cutoff in TEST_CUTOFFS){
    main(study_short="meta_bcac_ukb_ovr", test=F, by_chrom=T, z_cutoff=cutoff)
	}
  # FOCUS_RESULTS %>% walk(~main(.))
  # FOCUS_RESULTS %>% walk(~main(., breast_filt=T))
  # paste0("intrinsic_subtype_", 1:5) %>% walk(~main(., test=F, by_chrom=T))
  # main(study_short="meta_bcac_ukb_ovr", test=F, by_chrom=T)
}
