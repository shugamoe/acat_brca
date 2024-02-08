# https://lab-notes.hakyimlab.org/post/2022-06-10-qgt-training/
# Need "pvaluess", "MAF", "N"

get_study_imputed_gwas <- function(study_name, dir="../input/metaxcan_gwas_imputed/"){
  require(tidyverse)
	require(data.table)
	require(glue)

	return(fread(file.path(dir, glue("{study_name}.txt.gz"))))


}

get_study_n <- function(study_name){
	return(
				 switch(study_name,
								"AD_Bellenguez_GWAS_NG_2022" = 111326 + 677663,
								"AD_GWAS_only" = 21982 + 41944,
								"AD_GWAS_GWAX_meta" =  242341,
								"meta_analysis_BCAC_UKB_ovr_brca" = 4480 + 17588 + 7333 + 42892 + 9655 + 45494,
								"meta_analysis_BCAC_CIMBA_erneg_brca" = 7784 + 7782 + 1630 + 1712 + 4480 + 17588 + 7333 + 42892 + 9655 + 45494,
								"intrinsic_subtype_5" = 8602 + 91477,
								"intrinsic_subtype_4" = 2884 + 91477,
								"intrinsic_subtype_3" = 6350 + 91477,
								"intrinsic_subtype_2" = 6427 + 91477,
								"intrinsic_subtype_1" = 45253 + 91477,
								"BCAC_Overall_BreastCancer_EUR" = 122977 + 105974,
								"BCAC_ERPOS_BreastCancer_EUR" = 69501 + 105974 
				 )
				 )
}


get_study_s <- function(study_name){
	return(
				 switch(study_name,
								"AD_Bellenguez_GWAS_NG_2022" = 111326/(677663 + 111326),
								"AD_GWAS_only" = 21982/(41944+21928),
								"AD_GWAS_GWAX_meta" =  242341,
								"meta_analysis_BCAC_UKB_ovr_brca" = (4480 + 7333 + 9655) / (17588 + 42892 + 45494 + 4480 + 7333 + 9655),
								"meta_analysis_BCAC_CIMBA_erneg_brca" = (7784 + 1630 + 4480 + 7333 + 9655) / (7782 + 1712 + 17588 + 42892 + 45494 + 7784 + 1630 + 4480 + 7333 + 9655),
								"intrinsic_subtype_5" = 8602 / (91477 + 8602),
								"intrinsic_subtype_4" = 2884 / (91477 + 2884),
								"intrinsic_subtype_3" = 6350 / (91477 + 6350),
								"intrinsic_subtype_2" = 6427 / (91477 + 6427),
								"intrinsic_subtype_1" = 45253 / (91477 + 45253),
								"BCAC_Overall_BreastCancer_EUR" = 122977 / (105974 + 122977),
								"BCAC_ERPOS_BreastCancer_EUR" = 69501 / (105974 + 69501)
				 )
				 )
}

do_coloc_abf_for_chr_tiss <- function(study, tiss, gwas_df, cur_tiss_chr_parquet, coloc_do_by_df, eqtl_or_sqtl){
	require(coloc)
	require(glue)
	require(tidyverse)
	require(data.table)

	do_one_sqtl <- function(gene, intron){
	  coloc_sqtl <- cur_tiss_chr_parquet %>%
			filter(gene == !!gene,
						 intron == !!intron) %>%
			rename(pvalues = pval_nominal) %>%
			filter(pvalues > 0, N > 0) %>%
		  drop_na(N) %>%
		  as.list()

		gwas_s <- gwas_df %>%
			pull(s) %>% head(1)

		coloc_gwas <- gwas_df %>%
			filter(snp %in% coloc_sqtl$snp,
						 pvalues > 0) %>%
			select(-s)
		if (nrow(coloc_gwas) == 0){
			print(glue("{gene}:{intron} has no overlapping gwas SNPs"))
		  return(tibble())
		}
	  coloc_gwas <- coloc_gwas %>%
			as.list()

		coloc_sqtl[["type"]] <- "quant"
		coloc_gwas[["type"]] <- "cc"
		coloc_gwas[["s"]] <- gwas_s

		coloc_res <- runonce::save_run({coloc.abf(coloc_gwas, coloc_sqtl)},
			glue("../output/intermediate/{study}__{tiss}__{gene}__{intron}.rds"))

		best_h4_snp <- coloc_res$results %>% arrange(desc(`SNP.PP.H4`)) %>% head(1)

		return_tib <- tibble(gene = gene, intron = intron, tiss = tiss, H4 = coloc_res$summary[["PP.H4.abf"]],
												 top_h4_snp = best_h4_snp$snp, top_h4_snp_post=best_h4_snp$`SNP.PP.H4`)
		return(return_tib)
	}
	do_one_eqtl <- function(gene){
	  coloc_eqtl <- cur_tiss_chr_parquet %>%
			filter(gene == !!gene) %>%
			rename(pvalues = pval_nominal) %>%
			filter(pvalues > 0, N > 0) %>%
		  drop_na(N) %>%
		  as.list()

		gwas_s <- gwas_df %>%
			pull(s) %>% head(1)

		coloc_gwas <- gwas_df %>%
			filter(snp %in% coloc_eqtl$snp,
						 pvalues > 0) %>%
			select(-s)
		if (nrow(coloc_gwas) == 0){
			print(glue("{gene} has no overlapping gwas SNPs"))
		  return(tibble())
		}
	  coloc_gwas <- coloc_gwas %>%
			as.list()

		coloc_eqtl[["type"]] <- "quant"
		coloc_gwas[["type"]] <- "cc"
		coloc_gwas[["s"]] <- gwas_s

		coloc_res <- runonce::save_run({coloc.abf(coloc_gwas, coloc_eqtl)},
			glue("../output/intermediate/{study}__{tiss}__{gene}.rds"))

		best_h4_snp <- coloc_res$results %>% arrange(desc(`SNP.PP.H4`)) %>% head(1)

		return_tib <- tibble(gene = gene, tiss = tiss, H4 = coloc_res$summary[["PP.H4.abf"]],
												 top_h4_snp = best_h4_snp$snp, top_h4_snp_post=best_h4_snp$`SNP.PP.H4`)
		return(return_tib)
	}

	if (eqtl_or_sqtl == "eqtl"){
	  return_res <- coloc_do_by_df %>% pmap_dfr(do_one_eqtl)
	} else {
	  return_res <- coloc_do_by_df %>% pmap_dfr(do_one_sqtl)
	}
	return(return_res)
}

rollup_coloc_results <- function(results, study, eqtl_or_sqtl, 
																 tiss_lists=paste0("../input/tissue_lists/",
																									 c("11tiss.txt"))){
	genes_unique <- results %>%
		pull(gene) %>% unique()

	rollup_one_tisslist <- function(tiss_list_path){
		want_tissues <- fread(tiss_list_path, header=F)[["V1"]]
		results_filtered <- results %>%
			filter(tiss %in% want_tissues)
    tiss_list_basename_no_ext <- basename(tiss_list_path) %>%
      str_replace_all("\\.txt", "")

		results_filtered <- results %>%
			filter(tiss %in% want_tissues)

		get_max_tissue <- function(gname){
		  g_results <- results_filtered %>%
				filter(gene == !!gname)

			max_h4_df <- g_results %>%
				group_by(gene) %>%
				summarize(max_h4 = max(H4, na.rm=T))

			max_h4 <- max(max_h4_df$max_h4, na.rm=T)

			gmax_results <- g_results %>%
				filter(H4 == max_h4_df$max_h4)

			if (nrow(gmax_results) > 1){
				print(glue("{gname} has a repeat in {study} {eqtl_or_sqtl}, be sure to check it"))
				combined_tiss_names <- paste(unique(gmax_results$tiss) %>% sort(), collapse=",")
				if (eqtl_or_sqtl == "sqtl"){
				  combined_intron_names <- paste(unique(gmax_results$intron) %>% sort(), collapse=",")
				}
				combined_top_h4_snp <- paste(unique(gmax_results$top_h4_snp) %>% sort(), collapse=",")
				combined_top_h4_snp_post <- max(gmax_results$top_h4_snp_post)

				return_tib <- tibble(gene = gname,
														 coloc_top_h4 = max_h4,
														 coloc_top_h4_tissue = combined_tiss_names,
														 coloc_top_h4_snp = combined_top_h4_snp,
														 coloc_top_h4_snp_post = combined_top_h4_snp_post)
				if (eqtl_or_sqtl == "sqtl"){
          return_tib <- return_tib %>%
						mutate(coloc_top_h4_intron = combined_intron_names)
				}
			} else {
				return_tib <- tibble(gene = gname,
														 coloc_top_h4 = max_h4,
														 coloc_top_h4_tissue = gmax_results$tiss,
														 coloc_top_h4_snp = gmax_results$top_h4_snp,
														 coloc_top_h4_snp_post = gmax_results$top_h4_snp_post)
				if (eqtl_or_sqtl == "sqtl"){
				  return_tib <- return_tib %>%
						mutate(coloc_top_h4_intron = gmax_results$intron)
				}
			}
			return_tib
		}

		write_df <- tibble(gname = genes_unique) %>%
			pmap_dfr(get_max_tissue)
		write_tsv(write_df, file.path("../output/final/", glue("gene_rolled_{eqtl_or_sqtl}_{tiss_list_basename_no_ext}_{study}.tsv")))
	}
  tibble(tiss_list_path=tiss_lists) %>%
		pmap(~rollup_one_tisslist(.))
}



main <- function(want_study, want_gene_file = "../input/brca_want_genes.txt", eqtl_or_sqtl="eqtl",
								 want_tiss_file="../input/tissue_lists/11tiss.txt"){
	require(tidyverse)
	require(data.table)
	require(glue)
	want_genes <- fread(want_gene_file, header=F)[["V1"]]
	want_tissues <- fread(want_tiss_file, header=F)[["V1"]]
	study_n <- get_study_n(want_study)
	study_s <- get_study_s(want_study)

	if (eqtl_or_sqtl == "sqtl"){
	  qtl_tag <- "sqtl_"
	} else {
		qtl_tag <- ""
	}

	if (eqtl_or_sqtl == "sqtl"){
		gene_intron_key <- fread("../input/combine_phenotype_groups.txt.gz") %>%
			filter(gene_id %in% want_genes) %>%
			as_tibble()
		gene_intron_key[c("tissue", "intron")] <- stringr::str_split_fixed(gene_intron_key$intron_id, pattern="\\.", 2)

		gene_intron_key <- gene_intron_key %>%
			filter(tissue %in% want_tissues) %>%
			separate(gtex_intron_id, c(NA, "n1", "n2", NA,"gene"), ":", convert=TRUE)

		want_ids <- gene_intron_key %>%
			pull(intron) %>% unique()
		  
	} else {
	  want_ids <- want_genes
	}

  gwas_df <- get_study_imputed_gwas(want_study) %>%
		mutate(N = study_n,
					 s = study_s) %>%
		select(all_of(c("panel_variant_id", "pvalue", "N", "s", "frequency"))) %>%
		rename(snp = panel_variant_id,
					 MAF = frequency, 
					 pvalues=pvalue)

	res_tib <- tibble()
	for (cur_tiss in want_tissues){ # debug
	  did_eqtl_tiss <- F
	  for (cur_chr in 1:22){ # debug
			if (isTRUE(did_eqtl_tiss)){
			  next
			}
			print(glue("Tissue: {cur_tiss} | Chr {cur_chr} | {eqtl_or_sqtl}"))
			if (eqtl_or_sqtl == "eqtl"){
			  fname <- glue("{cur_tiss}.v8.allpairs.parquet")
			} else {
			  fname <- glue("{cur_tiss}.v8.cis_sqtl.all_pairs.chr{cur_chr}.parquet")
			}
			fpath <- file.path("../input/parquet_raw_qtl/", fname)
		  # cur_tiss_chr_parquet <- arrow::read_parquet(glue("../input/parquet_raw_qtl/{cur_tiss}.v8.EUR.{qtl_tag}allpairs.chr{cur_chr}.parquet"))
			if (eqtl_or_sqtl == "eqtl"){
				og_parquet <- arrow::read_parquet(fpath)
				cur_tiss_chr_parquet <- og_parquet %>%
					rename(gene = gene_id) %>%
					filter(gene %in% want_genes)
			} else {
		    cur_tiss_chr_parquet <- arrow::read_parquet(fpath)
			}
			if (eqtl_or_sqtl == "sqtl"){
				cur_tiss_chr_parquet <- cur_tiss_chr_parquet %>%
					filter(str_detect(phenotype_id, paste0(want_genes, collapse="|")))
			  if (nrow(cur_tiss_chr_parquet) == 0){
			  	next
			  }
				cur_tiss_chr_parquet <- cur_tiss_chr_parquet %>%
					separate(phenotype_id, c(NA, "n1", "n2", NA, "gene"), ":", convert=TRUE) %>%
					left_join(gene_intron_key %>% filter(tissue == !!cur_tiss) %>%
										select(n1, n2, gene, intron, tissue), by = c("n1",
																																			"n2",
																																			"gene")) %>%
					filter(gene %in% want_genes)

        coloc_do_by_df <- cur_tiss_chr_parquet %>%
					group_by(gene, intron) %>%
					summarize(n = n()) %>%
					select(-all_of(c("n")))
			} else {
				coloc_do_by_df <- cur_tiss_chr_parquet %>%
					group_by(gene) %>%
					summarize(n = n()) %>%
					select(-all_of(c("n")))
			}

			cur_tiss_chr_parquet <- cur_tiss_chr_parquet %>%
				mutate(N = as.integer((1 / af) * ma_samples)) %>%
				rename(snp = variant_id,
							 MAF=af)

			cur_tiss_chr_res <- do_coloc_abf_for_chr_tiss(want_study, cur_tiss, gwas_df, cur_tiss_chr_parquet, coloc_do_by_df, eqtl_or_sqtl)
			res_tib <- bind_rows(res_tib, cur_tiss_chr_res)
			if (eqtl_or_sqtl == "eqtl"){
			  did_eqtl_tiss <- T
			}
		}
	}
	print(glue("Rolling up coloc results for {want_study}"))
	write_tsv(res_tib, file.path("../output/final", glue("unrolled_{eqtl_or_sqtl}__{want_study}.tsv")))
	rolled_res_tib <- rollup_coloc_results(res_tib, want_study, eqtl_or_sqtl)
}



# STUDIES <- c("AD_Bellenguez_GWAS_NG_2022", "AD_GWAS_only", "AD_GWAS_GWAX_meta")
library(tidyverse)
STUDIES <- c("meta_analysis_BCAC_UKB_ovr_brca", "meta_analysis_BCAC_CIMBA_erneg_brca",
						 "intrinsic_subtype_5", "intrinsic_subtype_4", "intrinsic_subtype_3", "intrinsic_subtype_2",
						 "intrinsic_subtype_1", "BCAC_Overall_BreastCancer_EUR", "BCAC_ERPOS_BreastCancer_EUR")
EQTL_OR_SQTL <- c("eqtl", "sqtl")
param_grid <- expand.grid(want_study = STUDIES, eqtl_or_sqtl = c("eqtl", "sqtl"), stringsAsFactors=F)
if (interactive()){
	# main(STUDIES[1], eqtl_or_sqtl="sqtl")
	# main(STUDIES[1], eqtl_or_sqtl="eqtl")
} else {
	# main(STUDIES[1], eqtl_or_sqtl="sqtl")
	# main(STUDIES[1], eqtl_or_sqtl="eqtl")
	args = commandArgs(trailingOnly=TRUE)
	print(args)
	main(STUDIES[as.integer(args[1])], eqtl_or_sqtl=EQTL_OR_SQTL[as.integer(args[2])])
	# param_grid %>% pwalk(main)
}
