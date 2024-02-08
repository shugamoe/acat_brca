get_model_weight_source <- function(tissue, 
                                    eqtl_or_sqtl
                                    ){
	require(glue)
  db_dir <- glue("../input/gtex_v8_{eqtl_or_sqtl}_dbs_mashr/")
  db_weight_source <- file.path(db_dir, glue("mashr_{tissue}.db"))

  return(db_weight_source)
}

transform_sum_dat <- function(study, zscore_from_beta_se=F){
  require(data.table)
  require(tidyverse)
  require(glue)

  # study can be: schwartz_gwas, schwartz_gwasx, bellen, wight
  study_dir <- "../input/metaxcan_gwas_imputed"
  study_fpath <- switch(study, 
                        "bcac_erp" = file.path(study_dir, "BCAC_ERPOS_BreastCancer_EUR.txt.gz"),
                        "meta_bcac_cimba_erneg" = file.path(study_dir, "meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz"),
                        "bcac_ovr" = file.path(study_dir, "BCAC_Overall_BreastCancer_EUR.txt.gz"),
                        "meta_bcac_ukb_ovr" = file.path(study_dir, "meta_analysis_BCAC_UKB_ovr_brca.txt.gz"),
                        "intrinsic_subtype_1" = file.path(study_dir, "intrinsic_subtype_1.txt.gz"),
                        "intrinsic_subtype_2" = file.path(study_dir, "intrinsic_subtype_2.txt.gz"),
                        "intrinsic_subtype_3" = file.path(study_dir, "intrinsic_subtype_3.txt.gz"),
                        "intrinsic_subtype_4" = file.path(study_dir, "intrinsic_subtype_4.txt.gz"),
                        "intrinsic_subtype_5" = file.path(study_dir, "intrinsic_subtype_5.txt.gz")
                        )

  study_df <- fread(study_fpath,
                    select=c("panel_variant_id", # usually has rsid, except wightmen
                             "pvalue",
                             "chromosome", # of form chr<num>
                             "position",
                             "effect_allele",
                             "non_effect_allele",
                             "effect_size",
                             "standard_error",
                             "zscore")) %>%
    mutate(chromosome = as.integer(str_extract(chromosome, "\\d{1,2}"))) %>%
    rename(p_value = pvalue,
           base_pair_location = position,
           other_allele = non_effect_allele,
           beta = effect_size,
           SNP_ID = panel_variant_id
           )

  if (isTRUE(zscore_from_beta_se) & study != "wight"){
    study_df <- study_df %>%
      select(-all_of(c("zscore"))) %>%
      mutate(zscore = beta / standard_error)
  }

  final_cols <- c("SNP_ID", "p_value", "chromosome",  "base_pair_location",
                  "effect_allele", "other_allele", "beta", "standard_error", "zscore")

  # Wightmen had Inf/-Inf zscore, replace with approximately the most extreme zscore values seen
  study_df <- study_df %>%
    select(all_of(final_cols)) %>%
    mutate(zscore = ifelse(zscore == Inf, 40, zscore)) %>%
    mutate(zscore = ifelse(zscore == -Inf, -40, zscore)) %>%
    filter(!is.na(zscore)) # Will error out otherwise in impute_expr_z if na zscore

  return(study_df)
}
