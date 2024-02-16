GWAS_TYPES <- c("ERPOS", "ERNEG")
COND_TYPES <- c("eqtl", "sqtl", "both")
CLUMP_R2_VALS <- c(.65, .5, .4, .3, .2, .1)
TISS_PREFIXES <- c("11tiss")

make_results_long <- function(gwas_type="ERPOS"){
	require(tidyverse)
	require(data.table)
	require(glue)
	require(readxl)

	#if (file.exists(glue("../output/results_long/{gwas_type}_all_results.tsv"))){
		#return()
	# }

  gwas_type_master_df <- fread(glue("../output/{gwas_type}_master.tsv"))

	uniq_loci <- unique(gwas_type_master_df$Locus)

	final_gwas_type_long_df <- tibble()
	for (cloci in uniq_loci){
		for (tiss_prefix in TISS_PREFIXES){
			for (cond_type in COND_TYPES){
				for (r2_val in c(NA, CLUMP_R2_VALS)){
					base_cojo_run_id <- glue("{tiss_prefix}__{gwas_type}__{cloci}__gtex_{cond_type}")
					if (is.na(r2_val)){
						cojo_dir <- "../output/intermediate_data/raw_cojo"
						cojo_run_id <- base_cojo_run_id
					} else {
						cojo_dir <- "../output/intermediate_data/clumped_cojo"
						cojo_run_id <- glue("{base_cojo_run_id}__clump_r2_{r2_val}")
					}
					gwas_type_master_region <- gwas_type_master_df %>%
						filter(Locus == cloci, !is.na(`GWAS SNP`))
					results_path <- file.path(cojo_dir, glue("{cojo_run_id}.cma.cojo")) 
					olap_path <- file.path("../output/intermediate_data/overlap_tracking/", glue("{base_cojo_run_id}.tsv"))

					# If a SNP exists in the conditioning set and we want conditional effects for it, then it
					if (file.exists(results_path)){
						cojo_res_df <- fread(results_path) %>%
							mutate(z = b / se,
										 zC = bC / bC_se,
										 pC = 2 * pnorm(-abs(bC / bC_se))
										 ) %>%
							select(`SNP`, `b`, `bC`, `z`, `zC`, `p`, `pC`)
						if (file.exists(olap_path)){
							olap_df <- fread(olap_path)
							# The actual SNP list output from COJO that was used for conditioning
							actual_cond_df <-fread(file.path(cojo_dir, glue("{cojo_run_id}.given.cojo")))

							gwas_snps_in_cond_df <- actual_cond_df %>%
								filter(`SNP` %in% olap_df$chromosome_position) %>%
								select(`SNP`) %>%
								mutate(in_conditioning_set = TRUE)

							# Any overlap between the olap_df and actual_cond_df means that there
							# was a GWAS SNP that we want conditional effects for that actually
							# turned out to be in the set of predictive SNPs we want to condition
							# on. When this happens, we track the offending SNPs and be sure to
							# to assign an effect size. COJO should not be calculating
							# conditional effects for that SNP anyway though too.
							if (nrow(cojo_res_df) > 0 & nrow(gwas_snps_in_cond_df) > 0){
								add_to_final <- bind_rows(cojo_res_df, gwas_snps_in_cond_df) %>%
									mutate(`in_conditioning_set` = case_when(`in_conditioning_set` == F ~ NA,
																													 TRUE ~ `in_conditioning_set`))
							} else if (nrow(gwas_snps_in_cond_df) > 0){
								add_to_final <- gwas_snps_in_cond_df %>%
									mutate(`in_conditioning_set` = case_when(`in_conditioning_set` == F ~ NA,
																													 TRUE ~ `in_conditioning_set`))
							} else if (nrow(cojo_res_df) > 0){
								add_to_final <- cojo_res_df
							} else if (nrow(gwas_snps_in_cond_df) == 0 & nrow(cojo_res_df) == 0){
								next
							} else {
								message("What?")
								browser()
							}
						} else {
							# actual_cond_df <-fread(glue("../output/intermediate_data/raw_cojo/{cojo_run_id}.given.cojo"))
							# add_to_final <- bind_rows(cojo_res_df, actual_cond_df %>% select(`SNP`) %>% mutate(`in_geno_data` = TRUE))
							if (nrow(cojo_res_df) > 0){
								add_to_final <- cojo_res_df
							} else {
								next
							}

						}
						add_to_final <- add_to_final %>%
							mutate(tissue_set = tiss_prefix, 
										 clump_r2 = r2_val,
										 Locus = cloci,
										 in_geno_data = TRUE,
										 gtex_cond_type = cond_type
										 )

						# Check for mismatched by frequency SNPs
						bad_freq_path <- file.path(cojo_dir, glue("{cojo_run_id}.freq.badsnps"))
						if (file.exists(bad_freq_path)){
							bad_freq_df <- fread(bad_freq_path) %>%
								select(`SNP`) %>%
								mutate(freq_mismatch = T) %>%
								mutate(in_geno_data = TRUE)
							add_to_final <- add_to_final %>%
								bind_rows(bad_freq_df) %>%
								mutate(freq_mismatch = case_when(is.na(freq_mismatch) ~ F,
																								 freq_mismatch == TRUE ~ TRUE,
																								 TRUE ~ FALSE))

							num_mismatch <- sum(add_to_final$freq_mismatch)
							if (num_mismatch > 0){
							  message(glue("{cojo_run_id} has {num_mismatch} SNPs with mismatching frequencies"))
							}
						} else {
							add_to_final <- add_to_final %>%
								mutate(freq_mismatch=F)
						}

						if ("bC" %in% names(add_to_final)){
							if ("in_conditioning_set" %in% names(add_to_final)){
								add_to_final <- add_to_final %>%
										mutate(high_multivar_corr = case_when((is.na(`bC`) & !(freq_mismatch == TRUE) & !(in_conditioning_set == TRUE)) ~ TRUE,
																											TRUE ~ FALSE))
							} else {
								add_to_final <- add_to_final %>%
										mutate(high_multivar_corr = case_when((is.na(`bC`) & !(freq_mismatch == TRUE)) ~ TRUE,
																											TRUE ~ FALSE))
							}
							num_high_multivar_corr <- sum(add_to_final$high_multivar_corr)
							if (num_high_multivar_corr > 0){
								# message(glue("{cojo_run_id} has {num_high_multivar_corr} SNPs with high multivariate correlation with cond. SNPs"))
							}
						}

						add_to_final <- add_to_final %>%
							mutate(in_geno_data = case_when((in_geno_data == TRUE) ~ TRUE,
																							TRUE ~ FALSE))


						final_gwas_type_long_df <- bind_rows(final_gwas_type_long_df, add_to_final)
					}
				}
			}
		}
	}
	message(glue("Gathered results for {gwas_type}"))
	final_gwas_type_long_df <- final_gwas_type_long_df %>%
		mutate(high_multivar_corr = case_when(is.na(high_multivar_corr) ~ FALSE,
																					high_multivar_corr == TRUE ~ TRUE,
																					TRUE ~ FALSE),
					 in_conditioning_set = case_when(is.na(in_conditioning_set) ~ FALSE,
																					 in_conditioning_set == TRUE ~ TRUE,
																					 TRUE ~ FALSE)
					 )

	message(glue("{gwas_type}: {nrow(final_gwas_type_long_df)} rows before rsid join"))
	final_gwas_type_long_df <- final_gwas_type_long_df %>%
		left_join(gwas_type_master_df %>% select(`GWAS SNP`, `chromosome_position`) %>% distinct(),
							by = c("SNP" = "chromosome_position"))
	message(glue("{gwas_type}: {nrow(final_gwas_type_long_df)} rows after rsid join"))
  final_gwas_type_long_df <- final_gwas_type_long_df %>%
		select(`Locus`, `GWAS SNP`, `chromosome_position` = `SNP`, `tissue_set`,
					 `gtex_cond_type`, `clump_r2`, `b`, `bC`, `z`, `zC`, `p`, `pC`,
					 `in_geno_data`, `freq_mismatch`, `high_multivar_corr`,
					 `in_conditioning_set`) 
	if (nrow(final_gwas_type_long_df %>% filter(is.na(`GWAS SNP`))) > 0){
		final_gwas_type_long_df <- final_gwas_type_long_df %>%
			filter(!is.na(`GWAS SNP`) & !(freq_mismatch == TRUE))
		if (nrow(final_gwas_type_long_df %>% filter(is.na(`GWAS SNP`))) > 0){
			message("Rogue SNP, inspect `final_gwas_type_long_df`")
			browser()
		}
	}
	write_tsv(final_gwas_type_long_df, glue("../output/results_long/{gwas_type}_all_results.tsv"))
}

prelim_join_results <- function(gwas_type="ERPOS"){
	require(tidyverse)
	require(data.table)
	require(glue)

	master_res_long_df <- fread(glue("../output/results_long/{gwas_type}_all_results.tsv"))
  gwas_type_master_df <- fread(glue("../output/{gwas_type}_master.tsv"))

	return_df <- tibble()
	for (cloci in gwas_type_master_df$`Locus` %>% unique()){
    gwas_type_master_df <- fread(glue("../output/{gwas_type}_master.tsv"))
		master_res_long_region <- master_res_long_df %>%
			filter(Locus == cloci)

		if (nrow(master_res_long_region) == 0){
			message(glue("{gwas_type}: {cloci} has no COJO results"))
			next
    }

		# master_res_wide_region <- master_res_long_region %>%
		# 	pivot_wider(names_from=c(gtex_cond_type, clump_r2), values_from=c("z", "zC", "p", "pC", "in_geno_data", "in_conditioning_set"))
		# gwas_type_master_region <- gwas_type_master_df %>%
		# 	filter(Locus == cloci, !is.na(`GWAS SNP`))
		for (cond_type in COND_TYPES){
			region_type <- master_res_long_region %>%
				filter(gtex_cond_type == cond_type)
			if (nrow(region_type) == 0){
				message(glue("{gwas_type}: {cond_type}: {cloci} has no COJO results"))
				next
			}
			max_clump_r2 <- master_res_long_region %>%
				pull(clump_r2) %>%
				max()

			if (is.na(max_clump_r2)){
				region_type <- region_type %>%
					filter(is.na(clump_r2))
			} else {
				region_type <- region_type %>%
					filter(clump_r2 == max(clump_r2))
			}
			return_df <- bind_rows(return_df, region_type)
		}
	}

	write_tsv(return_df, glue("../output/results_long/{gwas_type}_no_clump_or_highest_clump_r2.tsv"))
	if (gwas_type == "ERPOS"){
	  og_master <- read_excel("../input/Additional_File_1_2024Feb13.xlsx", sheet="TABLE S3", range="A2:G791")
	} else if (gwas_type == "ERNEG"){
		og_master <- read_excel("../input/Additional_File_1_2024Feb13.xlsx", sheet="TABLE S4", range="A2:G251")
	}
	for (cond_type in COND_TYPES){
		# return_df_wide <- return_df %>%
		# 	pivot_wider(names_from=c(gtex_cond_type), values_from=c("zC", "pC", "in_conditioning_set")) %>%
		# 	mutate(`Locus_join_excel` = str_replace(Locus, "_", " "))
		# write_tsv(return_df_wide, glue("../output/results_wide/{gwas_type}_no_clump_or_highest_clump_r2.tsv"))
		return_df_ctype <- return_df %>%
			filter(gtex_cond_type == cond_type) %>%
			select(-gtex_cond_type)

		write_tsv(left_join(og_master, return_df_ctype %>% mutate(`Locus` = str_replace(Locus, "_", " "))), glue("../output/results_wide/joined_to_s3s4/{gwas_type}_gtex_{cond_type}_COJO_joined_to_original_sup_table.tsv"))
	}
}

if (interactive()){
	make_results_long("ERNEG")
	make_results_long()
	prelim_join_results("ERNEG")
	prelim_join_results()
} else {
	make_results_long("ERNEG")
	make_results_long()
	prelim_join_results("ERNEG")
	prelim_join_results()
}
