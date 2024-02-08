rollup_coloc_results <- function(results, study, eqtl_or_sqtl, 
																 tiss_lists=paste0("../input/tissue_lists/",
																									 c("11tiss.txt", "Breast_Mammary.txt"))){
	require(tidyverse)
	require(glue)
	require(data.table)

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
			# browser()
			return_tib
		}

		write_df <- tibble(gname = genes_unique) %>%
			pmap_dfr(get_max_tissue)
		write_tsv(write_df, file.path("../output/final/", glue("gene_rolled_{eqtl_or_sqtl}_{tiss_list_basename_no_ext}_{study}.tsv")))
	}
  tibble(tiss_list_path=tiss_lists) %>%
		pmap(~rollup_one_tisslist(.))
}

library(tidyverse)
library(glue)
STUDIES <- c("meta_analysis_BCAC_UKB_ovr_brca", "meta_analysis_BCAC_CIMBA_erneg_brca",
						 "intrinsic_subtype_5", "intrinsic_subtype_4", "intrinsic_subtype_3", "intrinsic_subtype_2",
						 "intrinsic_subtype_1", "BCAC_Overall_BreastCancer_EUR", "BCAC_ERPOS_BreastCancer_EUR")
EQTL_OR_SQTL <- c("eqtl", "sqtl")
if (interactive()){
	for (e_or_s in EQTL_OR_SQTL){
		for (cur_study in STUDIES){
			res_tib <- data.table::fread(glue("../output/final/unrolled_{e_or_s}__{cur_study}.tsv"))
			rollup_coloc_results(res_tib, cur_study, e_or_s)
		}
	}
} else {
}
