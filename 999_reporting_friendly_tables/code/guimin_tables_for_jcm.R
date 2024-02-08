# Makes the tables as specified in ../input/tablesjulian_format.docx
# Tons of pair-wise copy pasted (eqtl/sqtl) code.
source("dz_jl_tables_for_jcm.R") 
library(purrr)

GM_STUDIES <- c("meta_analysis_BCAC_UKB_ovr_brca")

make_guimin_tables <- function(
  study_name,
  tiss_filter = "../../02_acat_eqtl_sqtl/input/tissue_lists/11tiss.txt",
  gcode_convert = "../output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps_2023Jun02.tsv",
  guimin_table_names = c("guimin__{study_name}_table1_acat_splice_sig_1mb_from_index_snps.tsv",
                         "guimin__{study_name}_table2_acat_splice_sig_still_cojo_sig.tsv",
                         "guimin__{study_name}_supp_table_s1_splice_joint_or_breast_id.tsv",
                         "guimin__{study_name}_supp_table_s2__breast_id.tsv",
                         "guimin__{study_name}_table3_acat_splice_twas_spredixcan_drilldown.tsv"
                         ),
  odir="../output/guimin_tables/"
  ){
  require(tidyverse)

  table_df <- make_tables(study_name=study_name,
                         tiss_filter=tiss_filter,
                         gcode_convert=gcode_convert,
                         return_raw=T)
  table_df <- table_df %>%
    rename(max_pip_all_tissues = `sqtl.highest_pip`,
           tissue_with_max_pip = `sqtl.highest_pip_tissues`,
           pip_breast = `sqtl.breast_mammary_tissue_highest_pip`,
           max_pip_tissues_introns_in_cred_set = `sqtl.max_pip_num_in_cred_set`,
           breast_tissue_introns_in_cred_set = `sqtl.breast_mammary_tissue_introns_in_cred_set`,
           max_pip_all_tissues_ctwas = `ctwas.sqtl.highest_susie_pip`,
           tissue_with_max_pip_ctwas = `ctwas.sqtl.highest_susie_pip_tissues`,
           pip_breast_ctwas = `ctwas.sqtl.breast_mammary_tissue_highest_susie_pip`
           )

  # Use this table to add in a column "Reported in previous TWAS", to mark whether it was in our expr-TWAS
  expr_twas_sig <- table_df %>%
    filter(sig_eqtl_any == T)

  brca_known_vars_df <- fread("../input/Ten_TWAS_papers_gene_list_2023Feb28.tsv")
	browser()

	genes_to_add_to_ten_twas <- expr_twas_sig %>%
		rename(ENSG = `gene_id_no_decimal`,
					 pvalue_overall = `pval_eqtl`) %>%
		filter(!(`ENSG` %in% brca_known_vars_df$ENSG)) %>%
		mutate(`FirstAuthors` = 'Gao',
					 `PUBMEDID` = "1234567890") %>%
		select(all_of(c("ENSG", "FirstAuthors", "PUBMEDID", "pvalue_overall")))

	new_known_vars_df <- bind_rows(brca_known_vars_df, genes_to_add_to_ten_twas)
	write_tsv(new_known_vars_df, file.path("../output/", "Ten_TWAS_papers_and_meta_bcac_ukb_expr_twas_gene_list_2023Dec11.tsv"))
	browser()

	# Combine Ten TWAS with expr twas sig
  # missing <- c("ENSG00000116288.12", "ENSG00000116285.12")
  # ENSG00000116288.12
  # ENSG00000116285.12

  # Table 3 "eqtl version"
  # "Genes identified by splicing TWAS but not by expression"
  # List of columns (from Guimin):
  #  expression_ACAT_P           expression_most-siginicant-Tissue-z-score
  #  expression_most-siginicant-tissue_P    expression_most-siginicant-tissue
  #  expression_TWAS_breast_z-score           expression_TWAS_breast_P
  #
  #  most-significant-tissue_splicing_ACAT_P  most-significant-splicing-tissue
  #  introns            intron-s-Predixcan-z-score   intron-s-Predixcan-P
  gtable5_cols <- c("Gene Loci", "gene_name.v40", "genomic_coord", "pval_sqtl", "pval_sqtl_breast",
                    

                    # Below column block correspond to the column list above from Guimin
                    "pval_eqtl", "expression_most_sig_tissue_zscore",
                    "expression_most_sig_tissue_pval", "expression_most_sig_tissue",
                    "expression_TWAS_breast_zscore", "expression_TWAS_breast_P",
                    "splice_tissue_is_most_sig",
                    "tissue_splicing_ACAT_P", "most_significant_splicing_tissue_or_breast",
                    "intron", "intron_spredixcan_zscore", "intron_spredixcan_pval",

                    "gene_id.v26", "Distance from Index SNP(kb)", "thresh_sqtl", "thresh_sqtl_breast",
                    "Reported in our expression TWAS", "Reported in other TWAS")
  # Add eqtl columns for table3
  table5_genes <- table_df %>% filter(sig_sqtl_any == T, (is.na(sig_eqtl_any) | sig_eqtl_any == F)) %>% pull(gene_id.v26)
  zscore_eqtl_max_tissue_df <- get_zscore_detail(study_name, table5_genes,
    "eqtl", "tissue", tiss_filter
  ) %>%
    rename(`expression_most_sig_tissue_zscore` = zscore_for_max_z,
           `expression_most_sig_tissue_pval` = pval_for_max_z,
           `expression_most_sig_tissue` = eqtl_max_abs_zscore_tiss
            )

  zscore_eqtl_breast_df <- get_zscore_detail(study_name, table5_genes, 
    "eqtl", "table", tiss_filter
  ) %>% filter(tissue == "Breast_Mammary_Tissue") %>%
    rename(`expression_TWAS_breast_zscore` = zscore,
            `expression_TWAS_breast_P` = pvalue
            )

  gtable5_df <- table_df %>%
    filter(sig_sqtl_any == T, (is.na(sig_eqtl_any) | sig_eqtl_any == F)) %>%
    left_join(zscore_eqtl_max_tissue_df, by="gene_id.v26") %>%
    left_join(zscore_eqtl_breast_df, by="gene_id.v26") %>%
    mutate(`Reported in our expression TWAS` = ifelse(gene_id.v26 %in% expr_twas_sig$gene_id.v26, T, F)) %>%
    rename(`Reported in other TWAS` = brca_twas_reported) %>%
    mutate(`Reported in previous or our TWAS` = ifelse(
                   (`Reported in our expression TWAS` == T) | (`Reported in other TWAS` == T), T, F))
  print(glue("Guimin Tables: Size of table3 after adding eqtl columns: {nrow(gtable5_df)}"))

  # Add sqtl columns for table3

  # This function creates a table with ~ the columns below:
  #  most-significant-tissue_splicing_ACAT_P  most-significant-splicing-tissue
  #  introns            intron-s-Predixcan-z-score   intron-s-Predixcan-P
  create_splicing_breakdown_df <- function(){
    bdown_file_df <- fread(glue("../../02_acat_eqtl_sqtl/output/acat_tiss_breakdown_sqtl/11tiss__{study_name}_sqtl_acat_results.tsv")) %>%
      filter(group %in% table5_genes)

    get_single_gene_cols <- function(want_gene){
      # For splicing, "most significant tissue" means the tissue with the lowest *ACAT* p-value
      single_cols_want <- c("gene_id.v26", "splice_tissue_is_most_sig", "tissue_splicing_ACAT_P",
                            "intron", "most_significant_splicing_tissue_or_breast", "intron_spredixcan_zscore",
                            "intron_spredixcan_pval")
      filt_df <- bdown_file_df %>%
        filter(group == want_gene)
      filt_df <- filt_df %>% 
               filter((tissue_acat == min(tissue_acat, na.rm=T)) |
                       (tissue == "Breast_Mammary_Tissue")) %>%
        rename(`gene_id.v26` = group,
               `tissue_splicing_ACAT_P` = tissue_acat,
               `most_significant_splicing_tissue_or_breast` = tissue,
               `intron_spredixcan_zscore` = intron_zscore,
               `intron_spredixcan_pval` = intron_pval
               )
      min_acat <- filt_df %>% pull(`tissue_splicing_ACAT_P`) %>% min(., na.rm=T)
      filt_df <- filt_df %>%
        mutate(splice_tissue_is_most_sig = ifelse(`tissue_splicing_ACAT_P` == min_acat, TRUE, FALSE)) %>%
        select(all_of(single_cols_want))
      return(filt_df)
    }
    out_df <- bdown_file_df %>% pull(group) %>% unique() %>%
      map_dfr(~get_single_gene_cols(.))
    return(out_df)
  }
  splice_bdown_df <- create_splicing_breakdown_df()
  gtable5_df <- gtable5_df %>%
    left_join(splice_bdown_df, by="gene_id.v26") %>%
    select(all_of(gtable5_cols))
  print(glue("Guimin Tables: Size of table3 after adding sqtl columns: {nrow(gtable5_df)}"))
  write_delim(gtable5_df, file.path(odir, glue(guimin_table_names[5])), delim="\t")

  # Make/write Guimin tables

  # Table1 Splicing Significant and at least 1MB away from previous GWAS hits
  gtable1_cols <- c("Gene Loci", "gene_name.v40", "genomic_coord", "pval_sqtl", "pval_sqtl_breast", # "pip_focus_sqtl",

                    # These 3 columns are across tissue-exclusive models, a replacement for pip_focus_sqtl
                    "max_pip_all_tissues", "tissue_with_max_pip", "pip_breast",
                    "max_pip_all_tissues_ctwas", "tissue_with_max_pip_ctwas", "pip_breast_ctwas",

                    # "focus_introns_in_cred_set", 
                    "max_pip_tissues_introns_in_cred_set",
                    "breast_tissue_introns_in_cred_set",

                    
                    "enloc_sqtl_max_locus_rcp",
                    "enloc_sqtl_breast_max_locus_rcp",
                    "Reported in our expression TWAS", "Reported in other TWAS", 
                    "Reported in previous or our TWAS",
                    "gene_id.v26", "Distance from Index SNP(kb)", "thresh_sqtl", "thresh_sqtl_breast")
  gtable1_df <- table_df %>%
    filter(sig_sqtl_any == T,
           `Distance from Index SNP(kb)` > 1000) %>%
    mutate(`Reported in our expression TWAS` = ifelse(gene_id.v26 %in% expr_twas_sig$gene_id.v26, T, F)) %>%
    rename(`Reported in other TWAS` = brca_twas_reported) %>%
    mutate(`Reported in previous or our TWAS` = ifelse(
                   (`Reported in our expression TWAS` == T) | (`Reported in other TWAS` == T), T, F)) %>%
    select(any_of(gtable1_cols))
  write_delim(gtable1_df, file.path(odir, glue(guimin_table_names[1])), delim="\t")

  # Table2 Splicing Significant and at least 1MB away from previous GWAS hits
  gtable2_cols <- c("Gene Loci", "gene_name.v40", "genomic_coord", "pval_sqtl", "Index SNP", "Distance from Index SNP(kb)",
                    "pval_sqtl_cond", #"pip_focus_sqtl",

                    # These 3 columns are across tissue-exclusive models, a replacement for pip_focus_sqtl
                    "max_pip_all_tissues", "tissue_with_max_pip", "pip_breast",
                    "max_pip_all_tissues_ctwas", "tissue_with_max_pip_ctwas", "pip_breast_ctwas",

                    # "focus_introns_in_cred_set", 
                    "max_pip_tissues_introns_in_cred_set",
                    "breast_tissue_introns_in_cred_set",
                    
                    "enloc_sqtl_max_locus_rcp",
                    "Reported in our expression TWAS", "Reported in other TWAS",
                    "Reported in previous or our TWAS",
                    "gene_id.v26", "thresh_sqtl", "thresh_sqtl_breast")
  gtable2_df <- table_df %>%
                  filter(sig_sqtl_any == T,
                  sig_sqtl_any_cond == T) %>%
                  mutate(`Reported in our expression TWAS` = ifelse(gene_id.v26 %in% expr_twas_sig$gene_id.v26, T, F)) %>%
                  rename(`Reported in other TWAS` = brca_twas_reported) %>%
                  mutate(`Reported in previous or our TWAS` = ifelse(
                                  (`Reported in our expression TWAS` == T) | (`Reported in other TWAS` == T), T, F)) %>%
                  select(all_of(gtable2_cols))
  write_delim(gtable2_df, file.path(odir, glue(guimin_table_names[2])), delim="\t")

  # Supp Table1 Splicing Joint or Breast Significant
  gtable3_cols <- c("Gene Loci", "gene_name.v40", "gene_id.v26", "genomic_coord", "gene_type.v40", "pval_sqtl",
  "pval_sqtl_breast","pval_sqtl_cond", "pval_sqtl_breast_cond", "Distance from Index SNP(kb)",
  "Index SNP", "enloc_sqtl_max_locus_rcp", "enloc_sqtl_max_rcp_tissue", "enloc_sqtl_breast_max_locus_rcp",
  # "pip_focus_sqtl",

                    # These 3 columns are across tissue-exclusive models, a replacement for pip_focus_sqtl
                    "max_pip_all_tissues", "tissue_with_max_pip", "pip_breast",
                    "max_pip_all_tissues_ctwas", "tissue_with_max_pip_ctwas", "pip_breast_ctwas",
  
                    # "focus_introns_in_cred_set", 
                    "max_pip_tissues_introns_in_cred_set",
                    "breast_tissue_introns_in_cred_set",

  "Reported in our expression TWAS", "Reported in other TWAS", "Reported in previous or our TWAS",
  "FirstAuthors", "PUBMEDID", "sig_sqtl", "sig_sqtl_breast", "thresh_sqtl", "thresh_sqtl_breast"
  )
  brca_known_vars_df <- fread("../input/Ten_TWAS_papers_gene_list_2023Feb28.tsv") %>%
    select(all_of(c("ENSG", "FirstAuthors", "PUBMEDID")))
  gtable3_df <- table_df %>%
    filter((sig_sqtl == T) | (sig_sqtl_breast == T)) %>%
    mutate(`Reported in our expression TWAS` = ifelse(gene_id.v26 %in% expr_twas_sig$gene_id.v26, T, F)) %>%
    rename(`Reported in other TWAS` = brca_twas_reported) %>%
    mutate(`Reported in previous or our TWAS` = ifelse(
                    (`Reported in our expression TWAS` == T) | (`Reported in other TWAS` == T), T, F)) %>%
    left_join(brca_known_vars_df, by = c("gene_id_no_decimal" = "ENSG")) %>%
    select(all_of(gtable3_cols)) %>%
    arrange(sig_sqtl_breast, sig_sqtl) 
  write_delim(gtable3_df, file.path(odir, glue(guimin_table_names[3])), delim="\t")

  # Supp Table2 Breast Significant
  gtable4_cols <- c("gene_name.v40", "genomic_coord", "strand", "gene_id.v26", "gene_type.v40", "zscore_11tiss_abs_val_max",
  "pval_sqtl_breast", "pval_sqtl", 
  "Reported in our expression TWAS", "Reported in other TWAS", "Reported in previous or our TWAS",
  "FirstAuthors", "PUBMEDID", "sig_sqtl", "sig_sqtl_breast", "thresh_sqtl", "thresh_sqtl_breast"
  )
  # gene_id.v26`, `tissue`, `zscore`
  zscore_df <- get_zscore_detail(study_name, 
    table_df %>%
      filter(sig_sqtl_breast == T) %>% 
      pull("gene_id.v26"), "sqtl", "table", tiss_filter
  ) %>%
    mutate(abs_zscore = abs(zscore))
  zscore_df_grouped <- zscore_df %>%
   group_by(gene_id.v26) %>%
   summarize(abs_zscore_11tiss_abs_val_max = max(abs_zscore, na.rm=T))

  zscore_master <- left_join(zscore_df, zscore_df_grouped) %>%
    filter(abs_zscore == abs_zscore_11tiss_abs_val_max) %>%
    rename(zscore_11tiss_abs_val_max = zscore)

  gtable4_df <- table_df %>%
    left_join(zscore_master, by = c("gene_id.v26")) %>%
    filter(sig_sqtl_breast == T) %>%
    mutate(`Reported in our expression TWAS` = ifelse(gene_id.v26 %in% expr_twas_sig$gene_id.v26, T, F)) %>%
    rename(`Reported in other TWAS` = brca_twas_reported) %>%
    mutate(`Reported in previous or our TWAS` = ifelse(
                    (`Reported in our expression TWAS` == T) | (`Reported in other TWAS` == T), T, F)) %>%
    left_join(brca_known_vars_df, by = c("gene_id_no_decimal" = "ENSG")) %>%
    select(all_of(gtable4_cols)) %>%
    arrange(sig_sqtl_breast, sig_sqtl) %>%
    distinct(., gene_name.v40, genomic_coord, strand, gene_id.v26, gene_type.v40, pval_sqtl_breast, pval_sqtl, `Reported in our expression TWAS`, `Reported in other TWAS`, `Reported in previous or our TWAS`, FirstAuthors, PUBMEDID, sig_sqtl, sig_sqtl_breast, thresh_sqtl, thresh_sqtl_breast, .keep_all=TRUE) # Remove flipped values of `zscore_11tiss_abs_val_max`
  
  write_delim(gtable4_df, file.path(odir, glue(guimin_table_names[4])), delim="\t")

  return()
}

if (!interactive()) {
  make_guimin_tables(GM_STUDIES[1])
  # GM_STUDIES %>% walk(~make_tables(.,
  #                               tiss_filter = "../../02_acat_eqtl_sqtl/input/tissue_lists/11tiss.txt",
  #                               gcode_convert = "../output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps_2023Jun02.tsv",
  #                               table_names = paste0("from_dz_jl_", c("supp_table1__{study_name}.tsv",
  #                                               "table1_twas_unrep__{study_name}.tsv",
  #                                               "table2_cojo_sig__{study_name}.tsv",
  #                                               "table3_cand_causal__{study_name}.tsv",
  #                                               "supp_table2or3__{study_name}.tsv")),
  #                               odir="../output/guimin_tables/"
  #                               ))
} else {
  make_guimin_tables(GM_STUDIES[1])
  # GM_STUDIES %>% walk(~make_tables(.,
  #                               tiss_filter = "../../02_acat_eqtl_sqtl/input/tissue_lists/11tiss.txt",
  #                               gcode_convert = "../output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps.tsv",
  #                               table_names = paste0("from_dz_jl_", c("supp_table1__{study_name}.tsv",
  #                                               "table1_twas_unrep__{study_name}.tsv",
  #                                               "table2_cojo_sig__{study_name}.tsv",
  #                                               "table3_cand_causal__{study_name}.tsv",
  #                                               "supp_table2or3__{study_name}.tsv")),
  #                               odir="../output/guimin_tables/"
  #                               ))
}
