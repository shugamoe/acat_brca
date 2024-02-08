main <- function(){
	require(readxl)
	require(data.table)
	require(tidyverse)
	parent_gene_df <- tibble(parent_gene = c("FCGR1A", "FCGR2A", "FCGR3A",
																						"GBA", "GBA1", "F12", "PRSS1",
																						"TMPRSS2", "TMPRSS3", "ATP6AP1",
																						"PMS1", "PMS2", "RPL23A", "LCN1"))

	gcv40 <- read_csv("../input/gencode_v40_all.txt") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`gene_id`, "(ENSG\\d{1,})")) %>%
		select(`ens_gene_id_no_decimal`, `gene_name`)


	parent_gene_df <- parent_gene_df %>%
		left_join(gcv40, by = c("parent_gene" = "gene_name"))

	splice_df <- read_tsv("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/02_acat_eqtl_sqtl/output/acat_sqtl/11tiss__meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_results.tsv") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`group`, "(ENSG\\d{1,})")) %>%
		select(`ens_gene_id_no_decimal`, `multi_splice_twas_acat_pvalue` = `mtiss_acat`)

	splice_detail_df <- read_tsv("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/02_acat_eqtl_sqtl/output/acat_tiss_breakdown_sqtl/11tiss__meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_results.tsv") %>%
		filter(tissue == "Breast_Mammary_Tissue") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`group`, "(ENSG\\d{1,})")) %>%
		filter(ens_gene_id_no_decimal %in% parent_gene_df$`ens_gene_id_no_decimal`) %>%
		select(`ens_gene_id_no_decimal`,
					 `breast_tiss_splice_acat_pvalue` = `tissue_acat`) %>%
		distinct()

	express_acat_df <- read_tsv("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/02_acat_eqtl_sqtl/output/acat_eqtl/11tiss__meta_analysis_BCAC_UKB_ovr_brca_eqtl_acat_results.tsv") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`group`, "(ENSG\\d{1,})")) %>%
		filter(ens_gene_id_no_decimal %in% parent_gene_df$`ens_gene_id_no_decimal`) %>%
		select(`ens_gene_id_no_decimal`, `multi_express_twas_acat_pvalue` = `acat`)

	express_breast_df <- read_csv("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/01_spredixcan_eqtl_sqtl/output/spredixcan_eqtl_mashr/spredixcan_igwas_gtexmashrv8_meta_analysis_BCAC_UKB_ovr_brca__PM__Breast_Mammary_Tissue.csv") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`gene`, "(ENSG\\d{1,})")) %>%
		filter(ens_gene_id_no_decimal %in% parent_gene_df$`ens_gene_id_no_decimal`) %>%
		select(`ens_gene_id_no_decimal`, `breast_tiss_express_twas_pvalue` = `pvalue`,
					 `breast_tiss_express_twas_zscore` = `zscore`)

	parent_gene_df <- parent_gene_df %>%
		left_join(splice_df) %>%
		left_join(splice_detail_df) %>%
		left_join(express_acat_df) %>%
		left_join(express_breast_df)

	write_tsv(parent_gene_df, "../output/one_off/parent_gene_splice_and_express_info.tsv")
}


if (interactive()){
	main()
} else {
	main()
}
