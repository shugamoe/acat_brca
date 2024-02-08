# Combine CTWAS results by study, tissue, and whether eqtl or sqtl
arg_parser <- function(){
  require(argparse)
  parser <- ArgumentParser(description=paste0("",
                                              "",
                                              "",
                                              "",
                                              ""))

  parser$add_argument("--study", type="character", required=T, default=NULL,
                      help="Which study summary statistics to use."
  )
  parser$add_argument("--tissue_list", type="character", required=T, default="../input/tisuee_lists/11tiss.txt",
                      help="Which model (e.g. wingo_nc_2022) to use?"
  )
  parser$add_argument("--eqtl_or_sqtl", type="character", required=T, default=NULL,
                      help="Sourcing this model from the 'db' file or 'fusion' directly?"
  )
  parser_args  <- parser$parse_args()
  return(parser_args)
}

main <- function(study, tissue_list, eqtl_or_sqtl){
  require(data.table)
  require(tidyverse)
  require(glue)

  study_fpath <- switch(study, 
                        "bcac_erp" = file.path(study_dir, "BCAC_ERPOS_BreastCancer_EUR.txt.gz"),
                        "meta_bcac_cimba_erneg" = file.path(study_dir, "meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz"),
                        "bcac_ovr" = file.path(study_dir, "BCAC_Overall_BreastCancer_EUR.txt.gz"),
                        "meta_bcac_ukb_ovr" = file.path(study_dir, "meta_analysis_BCAC_UKB_ovr_brca.txt.gz"),
                        )
}

old_main <- function(){

  studies <- c("schwartz_gwas", "schwartz_gwasx", "bellen",
               "wight", "jansen")
  mod_types <- c("best", "lasso")

  col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0")

  col_types <- c("zscore", "pvalue", "pred_perf_r2", "bfr_thresh", "bfr_sig", "susie_pip", "mu2")
  for (study in studies){
    for (col_type in col_types){
      for (mod_type in mod_types){
        col_order <- append(col_order, glue("{study}.{mod_type}.{col_type}"))
      }
    }
  }


  base_df <- fread("../input/all_pwas_spxcan_lasso_vs_best_models.tsv") 
  for (study in studies){
    for (model_type in mod_types){
      mod_pip_col_name <- glue("{study}.{model_type}.susie_pip") 
      mod_mu2_col_name <- glue("{study}.{model_type}.mu2") 
      cur_df <- fread(glue("../output/ctwas_rss/{study}_wingo_nc_2022_fusion_{model_type}.susieIrss.txt")) %>%
        select(all_of(c("id", "type", "susie_pip", "mu2"))) %>%
        filter(type == "gene") %>%
        separate(id, c("ENSG ID", "Gene Symbol"), sep="\\.") %>%
        rename(!!mod_pip_col_name := `susie_pip`,
               !!mod_mu2_col_name := `mu2`) %>%
        select(-all_of(c("type")))

      base_df <- base_df %>%
        left_join(cur_df)
    }
  }

  base_df <- base_df %>%
    select(all_of(col_order)) %>%
    rowwise("ENSG ID") %>%
    mutate(best_sm_w_b_acat = single_acat(c(`schwartz_gwasx.best.pvalue`,
                                            `wight.best.pvalue`,
                                            `bellen.best.pvalue`)),
           best_sm_w_b_j_acat = single_acat(c(`schwartz_gwasx.best.pvalue`,
                                              `wight.best.pvalue`,
                                              `bellen.best.pvalue`,
                                              `jansen.best.pvalue`))
    ) %>% ungroup()

  write_delim(base_df, "../output/all_pwas_spxcan_and_ctwas_lasso_vs_best_models.tsv", delim="\t")
}

get_model_weight_source <- function(model, 
                                    model_source,
                                    method # Only accepts "lasso" or "best" now
                                    ){
  db_dir <- "../input/model_db_dir/"
  db_weight_source <- switch(model,
                             "wingo_nc_2022" = file.path(db_dir, glue("Wingo_NC_2022_pqtl_fusion_{method}_models_predixcan.db"))
                             )


  fusion_weight_source <- switch(model,
                                 "wingo_nc_2022" = "../input/Wingo_NC_2022/shared_mechanism.pQTLs.fusion.WEIGHTS_unzip/SampleID-map.nonNAfam.fusion.WEIGHTS/train_weights")

  if (model_source == "fusion"){
    return(fusion_weight_source)
  } else if (model_source == "db"){
    return(db_weight_source)
  }
}

# SNP_ID p_value chromosome base_pair_location effect_allele other_allele beta standard_error
# rs61769339 0.5322664744291312 1 662622 A G 0.012680341449 0.0203031660791
transform_sum_dat <- function(study, zscore_from_beta_se=T){
  require(data.table)
  require(tidyverse)
  require(glue)
  # final_cols <- 

  # study can be: schwartz_gwas, schwartz_gwasx, bellen, wight
  study_dir <- "../input"
  study_fpath <- switch(study, 
                        "schwartz_gwas" = file.path(study_dir, "AD_Schwartz_GWAS.txt.gz"),
                        "schwartz_gwasx" = file.path(study_dir, "AD_Schwartz_GWASX.txt.gz"),
                        "bellen" = file.path(study_dir, "AD_Bellenguez_GWAS_NG_2022.txt.gz"),
                        "wight" = file.path(study_dir, "AD_Wightmen_NG_2021.txt.gz"))

  if (study == "wight"){
    # 1	rs3094315	0.48877594	752566	G	A
    # 1	rs3131972	0.4888678	752721	A	G
    rsid_hg19_df <- fread("../input/master_chr37_rsid.bim",
                          header=F)
                          
    setnames(rsid_hg19_df, c("chromosome", "rsid", "cm", "base_pair_location", "a1", "a2"))

    rsid_hg19_df <- rsid_hg19_df %>%
      select(all_of(c("chromosome", "rsid", "base_pair_location")))
  }

  study_df <- fread(study_fpath,
                    select=c("variant_id", # usually has rsid, except wightmen
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
           SNP_ID = variant_id
           )

  if (isTRUE(zscore_from_beta_se) & study != "wight"){
    study_df <- study_df %>%
      select(-all_of(c("zscore"))) %>%
      mutate(zscore = beta / standard_error)
  }

  if (study == "wight"){
    study_df <- left_join(study_df, rsid_hg19_df)

    na_rsid_still <- nrow(study_df %>% filter(is.na(rsid)))

    study_df <- study_df %>%
      mutate(SNP_ID = rsid) %>%
      filter(!is.na(rsid))
  }

  final_cols <- c("SNP_ID", "p_value", "chromosome",  "base_pair_location",
                  "effect_allele", "other_allele", "beta", "standard_error", "zscore")

  # Wightmen had Inf/-Inf zscore, replace with approximately the most extreme zscore values seen
  study_df <- study_df %>%
    select(all_of(final_cols)) %>%
    mutate(zscore = ifelse(zscore == Inf, 40, zscore)) %>%
    mutate(zscore = ifelse(zscore == -Inf, -40, zscore))

  return(study_df)
}

single_acat <- function(pvals, weights=NULL, replace_pval = 1 - 1e-5){
  pvals[pvals == 1] <- replace_pval
  pvals <- pvals[!is.na(pvals)]

  if (is.null(weights)){
    weights <- 1 / length(pvals)
  }

  s = sum(weights * tan((0.5 - pvals) * pi))
  acat = 0.5 - atan(s / length(pvals)) / pi

  return(acat)
}

if (!interactive()){
  main()
} else {
  # main()
}
