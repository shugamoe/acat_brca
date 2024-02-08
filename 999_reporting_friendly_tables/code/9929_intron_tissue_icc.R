gene_icc <- function(){
	require(glue)
	require(tidyverse)
	require(data.table)
	require(psy)

	intron_df <- fread("../../02_acat_eqtl_sqtl/output/acat_tiss_breakdown_sqtl/11tiss__meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_results.tsv") %>%
		mutate(row=row_number()) %>%
		select(`gene_id.v26` = `group`, `intron`, `intron_zscore`, `tissue`) %>%
		distinct()

	model_df <- readRDS("/gpfs/data/gao-lab/Gao/manuscript_splicing_TWAS_acat/Haky files/gtex_extra_sqtl_from_bigquery_fewer_columns.RDS") %>%
		select(`intron` = `gene`, `tissue`, `pred_perf_R2`)

	gene_intron_df <- intron_df %>%
		left_join(model_df)
	write_tsv(gene_intron_df, "../output/one_off/11tiss_gene_intron_zscore_r2.tsv")
	return()

	genes <- c()
	consist_icc <- c(); agree_icc <- c(); r2_consist_icc <- c(); r2_agree_icc <- c(); z_p <- c(); r2_p <- c()
	unique_genes <- unique(gene_intron_df$`gene_id.v26`)
	for (gene in unique_genes){
	  cur_gene_df <- gene_intron_df %>%
			filter(`gene_id.v26` == gene)

		zscore_icc_input <- pivot_wider(cur_gene_df %>% select(-pred_perf_R2), names_from=tissue, values_from=intron_zscore, values_fn=max) %>%
			select(-all_of(c("gene_id.v26", "intron"))) %>%
			drop_na() %>%
			as.matrix()

		r2_icc_input <- pivot_wider(cur_gene_df %>% select(-intron_zscore), names_from=tissue, values_from=pred_perf_R2, values_fn=max) %>%
			select(-all_of(c("gene_id.v26", "intron"))) %>%
			drop_na() %>%
			as.matrix()

		zscore_icc_input_size = dim(zscore_icc_input)
		r2_icc_input_size = dim(r2_icc_input)

		z_p <- append(z_p, zscore_icc_input_size[2])
		r2_p <- append(r2_p, r2_icc_input_size[2])

		genes <- append(genes, gene)

		# N x p should be at least 2 x 2 (number of introns and number of tissue across)
		if (nrow(zscore_icc_input) >= 2 & dim(zscore_icc_input)[2] >= 2){
			zscore_res <- icc(zscore_icc_input)
			consist_icc <- append(consist_icc, zscore_res[[6]])
			agree_icc <- append(agree_icc, zscore_res[[7]])
		} else {
			consist_icc <- append(consist_icc,NA)
			agree_icc <- append(agree_icc, NA)
		}

		if (nrow(r2_icc_input) >= 2 & dim(r2_icc_input)[2] >= 2){
			r2_res <- icc(r2_icc_input)
			r2_consist_icc <- append(r2_consist_icc, r2_res[[6]])
			r2_agree_icc <- append(r2_agree_icc, r2_res[[7]])
		} else {
			r2_consist_icc <- append(r2_consist_icc,NA)
			r2_agree_icc <- append(r2_agree_icc, NA)
		} 
		if (length(genes) %% 100 == 0){
			message(glue("{length(genes)} / {length(unique_genes)}"))
			message(glue("zscore consist/agree {sum(!is.na(consist_icc))}"))
			message(glue("r2 consist/agree {sum(!is.na(r2_consist_icc))}"))
			message("")
		}
	}
	out_df <- tibble(genes=genes, zscore_consist_icc = consist_icc,
									 zscore_agree_icc=agree_icc,
									 r2_consist_icc=r2_consist_icc,
									 r2_agree_icc=r2_agree_icc,

									 # TODO(jcm) These 2 numbers are the same? Should they be?
									 zscore_p = z_p,
									 r2_p = r2_p)

	saveRDS(out_df, "../output/one_off/11tiss_zscore_and_r2_icc_by_gene_introns.RDS")

	# summary table
	# ugly gqq lol
	sum_df <- out_df  %>%
	 	mutate(across(contains("icc"), list(mean=~mean(.x, na.rm=T), sd=~ sd(.x, na.rm=T)), .names = "{.col}.{.fn}")) %>%
	 	select(-genes) %>%
	 	select(contains("icc.")) %>%
	 	distinct()
}

main <- function(){
	require(data.table)
	require(tidyverse)
	require(psy)
	require(glue)

	intron_df <- fread("../../02_acat_eqtl_sqtl/output/acat_tiss_breakdown_sqtl/11tiss__meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_results.tsv") %>%
		mutate(row=row_number()) %>%
		select(`intron`, `intron_zscore`, `tissue`) %>%
		distinct()

	model_r2_df <- fread("../output/one_off/intron_to_gene_pred_perf_R2_mapping.tsv") %>%
		select(`intron` = `intron_id`, `tissue`, `pred_perf_R2`, `gene_id`, `gene_name.v40`) %>%
		distinct()

	intron_r2_df <- intron_df %>%
		left_join(model_r2_df %>% select(-all_of(c("gene_id", "gene_name.v40"))))

	gene_intron_r2_df <- intron_df %>%
		left_join(model_r2_df)

	gene_icc(gene_intron_r2_df, "r2")
	gene_icc(gene_intron_r2_df, "zscore")

	browser()

	icc_input_df <- pivot_wider(intron_df, names_from=tissue, values_from=intron_zscore, values_fn = max) %>%
		select(-intron) %>%
		as.matrix()

	r2_icc_input_df <- pivot_wider(intron_r2_df %>% select(-intron_zscore), names_from=tissue, values_from=pred_perf_R2, values_fn = max) %>%
		select(-intron) %>%
		distinct() %>%
		as.matrix()

	# message(glue("Standard ICC: {date()}"))
	# icc_res <- icc(icc_input_df)

	# saveRDS(icc_res, "../output/one_off/11tiss_icc_intron.RDS")
	# message(glue("Saved to ../output/one_off/11tiss_icc_intron.RDS {date()}"))

	message(glue("Standard R2 ICC: {date()}"))
	r2_icc_res <- icc(r2_icc_input_df)

	saveRDS(r2_icc_res, "../output/one_off/11tiss_r2_icc_intron.RDS")
	message(glue("Saved to ../output/one_off/11tiss_r2_icc_intron.RDS {date()}"))

	# # Bootstrap
	# library(boot)
	# message(glue("Agreement bootstrap {date()}"))
  # icc.boot <- function(data,x) {icc(data[x,])[[7]]}
  # res <- boot(icc_input_df, icc.boot, 1000)
  # quantile(res$t,c(0.025,0.975)) # two-sided bootstrapped confidence interval of icc (agreement)
  # boot.ci(res,type="bca") # adjusted bootstrap percentile (BCa) confidence interval (better)

	# message(glue("Consistency bootstrap {date()}"))
  # icc.boot <- function(data,x) {icc(data[x,])[[6]]}
  # res <- boot(icc_input_df, icc.boot, 1000)
  # quantile(res$t,c(0.025,0.975)) # two-sided bootstrapped confidence interval of icc (agreement)
  # boot.ci(res,type="bca") # adjusted bootstrap percentile (BCa) confidence interval (better)

	# # Bootstrap
	# library(boot)
	# message(glue("R2 Agreement bootstrap {date()}"))
  # icc.boot <- function(data,x) {icc(data[x,])[[7]]}
  # res <- boot(r2_icc_input_df, icc.boot, 1000)
  # quantile(res$t,c(0.025,0.975)) # two-sided bootstrapped confidence interval of icc (agreement)
  # boot.ci(res,type="bca") # adjusted bootstrap percentile (BCa) confidence interval (better)

	# message(glue("R2 Consistency bootstrap {date()}"))
  # icc.boot <- function(data,x) {icc(data[x,])[[6]]}
  # res <- boot(r2_icc_input_df, icc.boot, 1000)
  # quantile(res$t,c(0.025,0.975)) # two-sided bootstrapped confidence interval of icc (agreement)
  # boot.ci(res,type="bca") # adjusted bootstrap percentile (BCa) confidence interval (better)
}


if (interactive()){
	# main()
	gene_icc()
} else {
	# main()
	gene_icc()
}
