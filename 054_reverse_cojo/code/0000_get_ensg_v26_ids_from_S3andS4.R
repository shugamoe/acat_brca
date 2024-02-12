main <- function(){
	require(tidyverse)
	require(data.table)
	require(readxl)
	require(glue)

	erpos_genes <- read_excel("../input/Supplementary_Tables_2023-07-02.xlsx", sheet="TABLE S3", range="A2:G791") %>%
		filter(`Region of TWAS or GWAS signals` == "Both") %>%
		filter(!is.na(`TWAS gene`)) %>%
		mutate(chr_num = str_extract(Locus, "^\\d{1,2}")) %>%
		mutate(chromosome = glue("chr{chr_num}")) %>%
		select(`chromosome`, `TWAS gene`) %>%
		distinct()

	erneg_genes <- read_excel("../input/Supplementary_Tables_2023-07-02.xlsx", sheet="TABLE S4", range="A2:G791") %>%
		filter(`Region of TWAS or GWAS signals` == "Both") %>%
		filter(!is.na(`TWAS gene`)) %>%
		mutate(chr_num = str_extract(Locus, "^\\d{1,2}")) %>%
		mutate(chromosome = glue("chr{chr_num}")) %>%
		select(`chromosome`, `TWAS gene`) %>%
		distinct()

	og_erpos_genes <- nrow(erpos_genes)
	og_erneg_genes <- nrow(erneg_genes)

	gname_ensg_df <- fread("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/999_reporting_friendly_tables/output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps.tsv") %>%
		select(`chromosome`, `gene_id.v26`, `gene_name.v40`, `gene_name.v38`, `gene_name.v26`) %>%
		distinct()

	erpos_genes_conv <- erpos_genes %>%
		left_join(gname_ensg_df %>% select(`chromosome`, `gene_name.v40`, `gene_id.v26`), by = c("TWAS gene" = "gene_name.v40", "chromosome")) %>% 
		filter(!is.na(`gene_id.v26`))

	erpos_genes_conv2 <- erpos_genes %>%
		left_join(gname_ensg_df %>% select(`chromosome`, `gene_name.v26`, `gene_id.v26`), by = c("TWAS gene" = "gene_name.v26", "chromosome")) %>% 
		filter(!is.na(`gene_id.v26`))

	erpos_genes_all_conv <- bind_rows(erpos_genes_conv, erpos_genes_conv2) %>%
		distinct()

	erpos_all_conv_num <- nrow(erpos_genes_all_conv)

	if (og_erpos_genes != erpos_all_conv_num){
		browser()
		error("Some erpos genes not converted")
	}

	erneg_genes_conv <- erneg_genes %>%
		left_join(gname_ensg_df %>% select(`chromosome`, `gene_name.v40`, `gene_id.v26`), by = c("TWAS gene" = "gene_name.v40", "chromosome")) %>% 
		filter(!is.na(`gene_id.v26`))

	erneg_genes_conv2 <- erneg_genes %>%
		left_join(gname_ensg_df %>% select(`chromosome`, `gene_name.v26`, `gene_id.v26`), by = c("TWAS gene" = "gene_name.v26", "chromosome")) %>% 
		filter(!is.na(`gene_id.v26`))

	erneg_genes_all_conv <- bind_rows(erneg_genes_conv, erneg_genes_conv2) %>%
		distinct()

	erneg_all_conv_num <- nrow(erneg_genes_all_conv)

	if (og_erneg_genes != erneg_all_conv_num){
		browser()
		error("Some erneg genes not converted")
	}

	write_csv(erpos_genes_all_conv, "../input/ERPOS_ensg_v26_ids.txt")
	write_csv(erneg_genes_all_conv, "../input/ERNEG_ensg_v26_ids.txt")
}
if (interactive()){
	main()
} else {
	main()
}
