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
  parser$add_argument("--tissue", type="character", required=T, default=NULL,
                      help="Which model (e.g. wingo_nc_2022) to use?"
  )
  parser$add_argument("--eqtl_or_sqtl", type="character", required=T, default=NULL,
                      help="Sourcing this model from the 'db' file or 'fusion' directly?"
  )

  parser$add_argument("--max_mag_snp_z", type="integer", required=T, default=20,
                      help="What is the highest allowable magnitude of SNP zscores?"
  )
  parser$add_argument("--filt_preharmo_z_snp", type="character", required=T, default=NULL,
                      help="RDS file for preharmnized (with LD) z_snp. zscore magnitude filtered."
  )
  parser$add_argument("--preharmo_z_snp", type="character", required=T, default=NULL,
                      help="RDS file for preharmnized (with LD) z_snp. zscore magnitude filtered."
  )

  parser$add_argument("--ncore", type="integer", required=T, default=NULL,
                      help="How many cores to use for par. execeution of ctwas."
  )
  parser$add_argument("--chr_num", type="integer", required=T, default=NULL,
                      help="What chromosome to run"
  )
  parser$add_argument("--thin", type="double", required=F, default=0.1,
                      help="thin argument"
  )
  
  parser_args  <- parser$parse_args()
  return(parser_args)
}


# Based off of /project2/guiming/xsun/1.brain_wingo/2.process_877_run_ctwas.R
# Modified with recent additions seen in /project2/guiming/xsun/1.brain_wingo_latest_sameld/ctwas.R
main <- function(study, tissue, eqtl_or_sqtl, chr_num, ncore, max_mag_snp_z, preharmo_z_snp, filt_preharmo_z_snp, thin){
  require(ctwas)
  require(tidyverse)
  require(glue)

  model_weight <- get_model_weight_source(tissue, eqtl_or_sqtl)

  message(glue("Using {model_weight} as model source."))

	# Get preharmonized/filtered z_snp
	z_snp_preharmo <- readRDS(preharmo_z_snp)$z_snp
	z_snp_filt_preharmo <- readRDS(filt_preharmo_z_snp)$z_snp


  # ld_R_dir <- paste0("../input/ctwas_gtex_v8_eur_ld")
  ld_R_dir <- paste0("../input/scratch_pdb_filt_ctwas_gtex_v8_eur_ld")


  z_res <- runonce::save_run({
                             impute_expr_z(z_snp=z_snp_preharmo, weight=model_weight, ld_R_dir = ld_R_dir,
                         method = method, outputdir="../output/impute_expr_z",
												 harmonize_z = F,
                         strand_ambig_action_z = "none", # Trust the predict.db models for strand ambig weights
                         strand_ambig_action_wgt = "none", 
                         outname = glue("{study}__max_mag_snp_z{max_mag_snp_z}__{tissue}__{eqtl_or_sqtl}_chr{chr_num}"),
                         ncore = ncore,
                         chrom = chr_num,
                         db_id="varID", # Hack to use the read the chr{num}_{pos}_{a0}_{a1}_b38 form of ID since that's what our ld ref uses
                         )},
                             file=glue("../output/impute_expr_z/{study}__max_mag_snp_z{max_mag_snp_z}__{tissue}__{eqtl_or_sqtl}_chr{chr_num}.rds"))
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

	ctwas_odir <- "../output/ctwas_rss/"
  ctwas_oname <- glue("{study}__max_mag_snp_z{max_mag_snp_z}__{tissue}__{eqtl_or_sqtl}_chr{chr_num}")
	param_calc_out <- tryCatch(
        {
            # Just to highlight: if you want to use more than one
            # R expression in the "try" part then you'll have to
            # use curly brackets.
            # 'tryCatch()' will return the last evaluated expression
            # in case the "try" part was completed successfully

            message("This is the 'try' part")
	          if (!file.exists(paste0(ctwas_odir, ctwas_oname, ".s2.susieIrssres.Rd"))){
              rv <- ctwas_rss(z_gene = z_gene,
						  				 	z_snp = z_snp_filt_preharmo, # Use preharmonized/snp z score magnitude filtered z_snp for parameter calculation
						  				 	ld_exprfs = ld_exprfs,
                  ld_regions_custom = glue("../input/regions/chr{chr_num}.txt"),
                  ld_R_dir = ld_R_dir,
                  thin = thin,
						  	 	max_snp_region = 20000,
						  	 	outputdir = ctwas_odir,
                  outname = ctwas_oname,
                  ncore = ncore,
						  	 	chrom=chr_num,
						  	 	merge=want_merge,
						  		fine_map=F
						  		)
						} else {
							rv <- "this estimation already finished"
						}
						rv
        },
        error=function(cond) {
            message("ERROR: Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            return(NULL)
        }
  )

	if (file.exists(paste0(ctwas_odir, ctwas_oname, ".s2.susieIrssres.Rd"))){
		print("Loading parameter estimation")
		load(paste0(ctwas_odir, ctwas_oname, ".s2.susieIrssres.Rd"))

		group_prior_rec <- group_prior_rec[,ncol(group_prior_rec)]
    group_prior_rec["SNP"] <- group_prior_rec["SNP"]*thin #adjust for thin argument

    group_prior_var_rec <- group_prior_var_rec[,ncol(group_prior_var_rec)]
    ctwas_rss(z_gene = z_gene, z_snp = z_snp_preharmo, ld_exprfs = ld_exprfs,
                ld_regions_custom = glue("../input/regions/chr{chr_num}.txt"),
                ld_R_dir = ld_R_dir,
                thin = thin,
	  					 	max_snp_region = 20000,
	  					 	outputdir = ctwas_odir,

								# Use parameters calculated with SNP set that was zscore
								# magnitude filtered
							  group_prior = group_prior_rec,
								group_prior_var = group_prior_var_rec,
								estimate_group_prior = F,
								estimate_group_prior_var = F,
								 
                outname = ctwas_oname,
                ncore = ncore,
								ncore_LDR = round(ncore/2),
	  					 	chrom=chr_num,
	  					 	merge=want_merge,
	  						fine_map=T)
		message("{ctwas_oname} completed")
	}
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
