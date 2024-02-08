main <- function(){
	require(readxl)
	require(data.table)
	require(tidyverse)
	table_s5 <- read_excel("../input/Supplemental_Tables__2023March10-Final-AJHG-expr.xlsx", sheet="TableS5", range="A2:J370")

	gcv26 <- read_tsv("../input/gencode_v26_all.txt") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`gene_id`, "(ENSG\\d{1,})")) %>%
		select(`ens_gene_id_no_decimal`,
					 `strand`)

	table_s2 <- read_excel("../input/Supplemental_Tables__2023March10-Final-AJHG-expr.xlsx", sheet="TableS2", range="A2:T311") %>%
    mutate(ens_gene_id_no_decimal = str_extract(`Ensembl Gene ID`, "(ENSG\\d{1,})")) %>%
		left_join(gcv26) %>%
		mutate(`Overall breast cancer risk` = "Yes",
					 `ER-` = NA,
					 `ER+` = NA,
					 `First Authors` = "Gao G",
					 `PUBMED ID` = "37164006",
					 `First Publication Date` = "6/1/2023"
					 ) %>%
		select(`Gene Symbol` = `Gene symbol`,
					 `hg38 position` = `Chromosome:Position (hg38)`,
					 `Gene type`,
					 `strand`,
					 `Overall breast cancer risk`,
					 `ER+`,
					 `ER-`,
					 `First Authors`,
					 `PUBMED ID`,
					 `First Publication Date`
					 )

	new_s5 <- table_s2 %>%
		filter(!(`Gene Symbol` %in% table_s5$`Gene Symbol`)) %>%
		bind_rows(table_s5, .)
		
	browser()

	write_tsv(new_s5, "../output/one_off/AJHG_s5_with_s2_non_overlapping_added.tsv")
}


if (interactive()){
	main()
} else {
	main()
}
