source("000_helpers.R")
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
  parser$add_argument("--max_mag_snp_z", type="integer", required=T, default=NULL,
                      help="What is the highest allowable magnitude of SNP zscores?"
  )
  parser_args  <- parser$parse_args()
  return(parser_args)
}


# Based off of /project2/guiming/xsun/1.brain_wingo/2.process_877_run_ctwas.R
# Modified with recent additions seen in /project2/guiming/xsun/1.brain_wingo_latest_sameld/ctwas.R
main <- function(study, max_mag_snp_z){
  require(ctwas)
  require(tidyverse)
  require(glue)

  study_df <- transform_sum_dat(study)
  # print(glue("zscore range: {range(study_df$zscore)}"))
  # return()

  # ld_R_dir <- paste0("../input/ctwas_gtex_v8_eur_ld")
  ld_R_dir <- paste0("../input/scratch_pdb_filt_ctwas_gtex_v8_eur_ld")

  # Do impute_expr_z 
  z_snp <- as.data.frame(cbind(study_df$SNP_ID, study_df$effect_allele, study_df$other_allele, study_df$zscore))
  colnames(z_snp) <- c("id", "A1", "A2", "z")
  z_snp$z <- as.numeric(z_snp$z)

	og_z_snp_num <- nrow(z_snp)

	z_snp <- z_snp %>%
		filter(abs(z) < max_mag_snp_z)

	filt_z_snp_num <- nrow(z_snp)

	print(glue("{og_z_snp_num - filt_z_snp_num} SNPs removed from {study} with abs(z) >= {max_mag_snp_z}"))
  # z_snp <- z_snp[!(z_snp$id %in% z_snp$id[duplicated(z_snp$id)]),] #drop multiallelic variants (id not unique)
  # ^^^New from: /project2/guiming/xsun/1.brain_wingo_latest_sameld/ctwas.R
  # For GTEx V8 models we can keep multiallelic variants since they do have
  # unique ids (not using rsid)
	preharmo_z <- runonce::save_run({preharmonize_z_ld(
    z_snp = z_snp,
    ld_R_dir = ld_R_dir,
    outputdir = "../output/preharmo/",
    outname = glue("{study}_max_mag_snp_z{max_mag_snp_z}"),
    harmonize_z = T,
    strand_ambig_action_z = "none",
    drop_multiallelic = F
    )}, file=glue("../output/preharmo/{study}_max_mag_snp_z{max_mag_snp_z}_save_run.rds")
	)
  warnings()
}

if (!interactive()){
  options(error = function() traceback(2))
  sink(stdout(), type="message")
  cli_args <- arg_parser()
  print(cli_args)
  do.call(main, cli_args)
} else {
}
