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
  parser$add_argument("--tissue", type="character", required=T, default=NULL,
                      help="Which model (e.g. wingo_nc_2022) to use?"
  )
  parser$add_argument("--eqtl_or_sqtl", type="character", required=T, default=NULL,
                      help="Sourcing this model from the 'db' file or 'fusion' directly?"
  )

  parser$add_argument("--ncore", type="integer", required=T, default=NULL,
                      help="How many cores to use for par. execeution of ctwas."
  )
  parser$add_argument("--chr_num", type="integer", required=T, default=NULL,
                      help="What chromosome to run"
  )
  
  parser_args  <- parser$parse_args()
  return(parser_args)
}


# Based off of /project2/guiming/xsun/1.brain_wingo/2.process_877_run_ctwas.R
# Modified with recent additions seen in /project2/guiming/xsun/1.brain_wingo_latest_sameld/ctwas.R
main <- function(study, tissue, eqtl_or_sqtl, chr_num, ncore){
  require(ctwas)
  require(tidyverse)
  require(glue)

  model_weight <- get_model_weight_source(tissue, eqtl_or_sqtl)

  message(glue("Using {model_weight} as model source."))

  study_df <- transform_sum_dat(study)
  # print(glue("zscore range: {range(study_df$zscore)}"))
  # return()

  # ld_R_dir <- paste0("../input/ctwas_gtex_v8_eur_ld")
  ld_R_dir <- paste0("../input/scratch_pdb_filt_ctwas_gtex_v8_eur_ld")

  # Do impute_expr_z 
  z_snp <- as.data.frame(cbind(study_df$SNP_ID, study_df$effect_allele, study_df$other_allele, study_df$zscore))
  colnames(z_snp) <- c("id", "A1", "A2", "z")
  z_snp$z <- as.numeric(z_snp$z)

  # z_snp <- z_snp[!(z_snp$id %in% z_snp$id[duplicated(z_snp$id)]),] #drop multiallelic variants (id not unique)
  # ^^^New from: /project2/guiming/xsun/1.brain_wingo_latest_sameld/ctwas.R
  # For GTEx V8 models we can keep multiallelic variants since they do have
  # unique ids (not using rsid)

  z_res <- runonce::save_run({
                             impute_expr_z(z_snp=z_snp, weight=model_weight, ld_R_dir = ld_R_dir,
                         method = method, outputdir="../output/impute_expr_z",
                         strand_ambig_action_z = "none", # Trust the predict.db models for strand ambig weights
                         recover_strand_ambig = F, # Don't need to flip weights based on correlations, should already match
                         outname = glue("{study}__{tissue}__{eqtl_or_sqtl}_chr{chr_num}"),
                         ncore = ncore,
                         chrom = chr_num,
                         db_id="varID", # Hack to use the read the chr{num}_{pos}_{a0}_{a1}_b38 form of ID since that's what our ld ref uses
                         wgt_keep_ambig=T # Hack to keep all weights from predict db file
                         )},
                             file=glue("../output/impute_expr_z/{study}__{tissue}__{eqtl_or_sqtl}_chr{chr_num}.rds"))
  print(glue("Imputation of zscores done for {tissue} {eqtl_or_sqtl} chr{chr_num}."))

  z_gene <- z_res$z_gene
  ld_exprfs <- z_res$ld_exprfs
  ld_exprfs <- ld_exprfs[!is.na(ld_exprfs)] # Remove NA stuff if doing one chrom at a time

  # Do ctwas_rss
  # Conditions to not merge
  if (tissue == "Vagina" & chr_num == 6 & eqtl_or_sqtl == "sqtl"){
    print("Not merging")
    want_merge <- F
  } else if (tissue == "Adipose_Subcutaneous" & chr_num == 6 & eqtl_or_sqtl == "eqtl"){
    print("Not merging")
    want_merge <- F
  } else {
    print("Merging regions.")
    want_merge <- T
  }

  ctwas_rss(z_gene = z_gene, z_snp = z_snp, ld_exprfs = ld_exprfs,
              ld_regions_custom = glue("../input/regions/chr{chr_num}.txt"),
              ld_R_dir = ld_R_dir,
              thin = 0.1, max_snp_region = 20000, outputdir = "../output/ctwas_rss",
              outname = glue("{study}__{tissue}__{eqtl_or_sqtl}_chr{chr_num}"),
              ncore = ncore, chrom=chr_num, merge=want_merge)
  warnings()
}

get_model_weight_source <- function(tissue, 
                                    eqtl_or_sqtl
                                    ){
  db_dir <- glue("../input/gtex_v8_{eqtl_or_sqtl}_dbs_mashr/")
  db_weight_source <- file.path(db_dir, glue("mashr_{tissue}.db"))

  return(db_weight_source)
}

# SNP_ID p_value chromosome base_pair_location effect_allele other_allele beta standard_error
# rs61769339 0.5322664744291312 1 662622 A G 0.012680341449 0.0203031660791
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

if (!interactive()){
  options(error = function() traceback(2))
  sink(stdout(), type="message")
  cli_args <- arg_parser()
  print(cli_args)
  do.call(main, cli_args)
} else {
}
