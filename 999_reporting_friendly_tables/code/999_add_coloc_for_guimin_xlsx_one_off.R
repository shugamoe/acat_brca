main <- function(){
	require(readxl)
	require(data.table)
	require(tidyverse)
	old_table <- read_excel("../input/Supp_table_10-7-23.xlsx", sheet="Table S1", range="A2:R252") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`Ensembl Gene ID`, "(ENSG\\d{1,})"))

	prev_twas_df <- fread("../input/Ten_TWAS_papers_gene_list_2023Feb28.tsv") %>%
		select(all_of(c("ENSG", "FirstAuthors", "PUBMEDID"))) %>%
		rename(ens_gene_id_no_decimal = ENSG)

	prev_expr_twas <- read_excel("../input/Supplemental_Tables__2023March10-Final.xlsx", sheet="TableS2",
															 range="A2:T311") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`Ensembl Gene ID`, "(ENSG\\d{1,})")) %>%
		mutate(in_table_s2 = TRUE) %>%
		select(`ens_gene_id_no_decimal`, `in_table_s2`)
	
	all_tiss_coloc <- fread("../../09_coloc/output/final/gene_rolled_sqtl_11tiss_meta_analysis_BCAC_UKB_ovr_brca.tsv") %>%
		select(`Ensembl Gene ID` = gene,
					 `sqtl.coloc_top_h4` = coloc_top_h4,
					 `sqtl.coloc_top_h4_tissue` = coloc_top_h4_tissue
					 )
	breast_coloc <- fread("../../09_coloc/output/final/gene_rolled_sqtl_Breast_Mammary_meta_analysis_BCAC_UKB_ovr_brca.tsv") %>%
		select(`Ensembl Gene ID` = gene,
					 `sqtl.coloc_top_h4_beast` = coloc_top_h4)

	# Get FOCUS in cred set
  focus_sqtl <- fread("../../999_reporting_friendly_tables/output/gene_rolled_focus_sqtl_11tiss_meta_analysis_BCAC_UKB_ovr_brca.focus.tsv") %>%
    select(all_of(c("ens_gene_id", "highest_pip", "introns_in_cred_set"))) %>%
    mutate(ens_gene_id_no_decimal = str_extract(`ens_gene_id`, "(ENSG\\d{1,})")) %>%
    rename(gene_id.v26 = ens_gene_id,
           pip_focus_sqtl = highest_pip,
           focus_introns_in_cred_set = introns_in_cred_set
           ) %>%
		select(`focus_introns_in_cred_set`, `ens_gene_id_no_decimal`)
	browser()
	 	
	new_table <- old_table %>%
		left_join(prev_expr_twas, by = "ens_gene_id_no_decimal") %>%
		left_join(focus_sqtl) %>%
		left_join(prev_twas_df) %>%
		left_join(all_tiss_coloc) %>%
		left_join(breast_coloc) %>%
		mutate(PUBMEDID = ifelse(in_table_s2 == TRUE & is.na(PUBMEDID), 37164006, PUBMEDID),
					 `FirstAuthors` = ifelse(in_table_s2 == TRUE & is.na(FirstAuthors), "Gao et al.", FirstAuthors)
					 ) %>%
		select(-all_of(c("ens_gene_id_no_decimal")))

	browser()

	write_tsv(new_table, "../output/one_off/Supp_table_10-23-23.tsv")
}


if (interactive()){
	main()
} else {
	main()
}
