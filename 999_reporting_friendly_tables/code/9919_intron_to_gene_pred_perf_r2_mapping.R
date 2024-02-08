main <- function(){
	require(data.table)
	require(tidyverse)

	model_df <- readRDS("/gpfs/data/gao-lab/Gao/manuscript_splicing_TWAS_acat/Haky files/gtex_extra_sqtl_from_bigquery_fewer_columns.RDS") %>%
		select(`intron_id` = `gene`, `tissue`, `pred_perf_R2`)

	map_df <- fread("../input/combine_phenotype_groups.txt.gz") %>%
		select(`intron_id`, `gene_id`, `cluster_id`) %>%
		separate(`intron_id`, sep="\\.", into = c("tissue", "intron_id"))

	gene_name_df <- fread("../output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps_2023Jun02.tsv") %>%
		select(`gene_id` = `gene_id.v26`, `gene_name.v40`)

	out_df <- left_join(model_df, map_df) %>%
		left_join(gene_name_df)

	out_df <- out_df %>%
		select(`gene_id`, `gene_name.v40`, `intron_id`, `tissue`, `cluster_id`, `pred_perf_R2`) %>%
		distinct()
	write_tsv(out_df, "../output/one_off/intron_to_gene_pred_perf_R2_mapping.tsv")
}


if (interactive()){
	main()
} else {
	main()
}
