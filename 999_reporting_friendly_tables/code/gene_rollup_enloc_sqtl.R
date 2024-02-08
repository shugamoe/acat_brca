ENLOC_RESULTS <- list.files(path="../../03_enloc/output/final_enloc_result_sqtl", pattern="*.txt",
                            full.names=TRUE)

library(tidyverse)
library(data.table)
library(glue)
GROUPINGS <- fread("../input/combine_phenotype_groups.txt.gz") %>% 
  as_tibble() %>% 
  rename(group = gene_id)
GROUPINGS[c("tissue", "intron")] <- stringr::str_split_fixed(GROUPINGS$intron_id, pattern="\\.", 2)

main <- function(input_enloc_sqtl_path,
                 tissue_filter="../../02_acat_eqtl_sqtl/input/tissue_lists/11tiss.txt",
                 odir="../output"){
  require(tidyverse)
  require(data.table)
  require(glue)
  require(purrr)

  print(input_enloc_sqtl_path)
  fname_only <- basename(input_enloc_sqtl_path)

  enloc_sqtl_df <- fread(input_enloc_sqtl_path) %>%
    rename(intron = molecular_qtl_trait) %>%
    as.tibble()

  enloc_sqtl_df <- left_join(GROUPINGS, enloc_sqtl_df, by = c("intron", "tissue")) %>%
    select(-intron_id)

  
  get_max_tissue_sqtl <- function(gname, enloc_sqtl, tissue_filter=NA){
    gene_df <- enloc_sqtl %>% filter(group == gname)

    # Get max breast tissue locus_rcp and corresponding intron(s)
    breast_df <- gene_df %>% filter(tissue == "Breast_Mammary_Tissue") %>%
      filter(locus_rcp == max(locus_rcp, na.rm=T)) %>%
      select(group, intron, locus_rcp, lead_coloc_SNP, lead_snp_rcp)
    if (nrow(breast_df) > 1){
      combined_intron_names <- paste(unique(breast_df$intron) %>% sort(), collapse=", ")
      combined_lead_coloc_SNP_names <- paste(unique(breast_df$lead_coloc_SNP) %>% sort(), collapse=", ")
      combined_lead_snp_rcp_names <- paste(unique(breast_df$lead_snp_rcp) %>% sort(), collapse=", ")
      breast_final <- tibble(group = gname, 
                             enloc_sqtl_breast_max_locus_rcp = breast_df$locus_rcp[1],
                             enloc_sqtl_breast_max_rcp_intron = combined_intron_names,
                             enloc_sqtl_breast_max_lead_coloc_SNP = combined_lead_coloc_SNP_names,
                             enloc_sqtl_breast_max_lead_snp_rcp = combined_lead_snp_rcp_names,
                             )
    } else if (nrow(breast_df) == 0){
      breast_final <- tibble(group = gname, 
                             enloc_sqtl_breast_max_locus_rcp = NA,
                             enloc_sqtl_breast_max_rcp_intron = NA,
                             enloc_sqtl_breast_max_lead_coloc_SNP = NA,
                             enloc_sqtl_breast_max_lead_snp_rcp = NA,
                             )
    } else {
      breast_df$lead_snp_rcp <- as.character(breast_df$lead_snp_rcp)
      breast_final <- breast_df %>% rename(enloc_sqtl_breast_max_rcp_intron = intron,
                                           enloc_sqtl_breast_max_locus_rcp = locus_rcp,
                                           enloc_sqtl_breast_max_lead_coloc_SNP = lead_coloc_SNP,
                                           enloc_sqtl_breast_max_lead_snp_rcp = lead_snp_rcp)
    }

    if (!is.na(tissue_filter)){
      tiss_filter <- fread(tissue_filter, header=F)
      gene_df <- gene_df %>% 
        filter(tissue %in% tiss_filter$V1)
    }

    # Get best all tissue introns, tissues, and locus rcps
    best_df <- gene_df %>%
      filter(locus_rcp == max(locus_rcp, na.rm=T)) %>%
      select(intron, locus_rcp, tissue, lead_coloc_SNP, lead_snp_rcp)

    min_rcp <- gene_df %>% pull(locus_rcp) %>% min(., na.rm=T)

    if (nrow(best_df) > 1){
      combined_tiss_names <- paste(unique(best_df$tissue) %>% sort(), collapse=", ")
      combined_intron_names <- paste(unique(best_df$intron) %>% sort(), collapse=", ")
      combined_lead_coloc_SNP_names <- paste(unique(best_df$lead_coloc_SNP) %>% sort(), collapse=", ")
      combined_lead_snp_rcp_names <- paste(unique(best_df$lead_snp_rcp) %>% sort(), collapse=", ")
      best_final <- tibble(enloc_sqtl_max_rcp_tissue = combined_tiss_names,
                           enloc_sqtl_max_locus_rcp = best_df$locus_rcp[1],
                           enloc_sqtl_min_locus_rcp = min_rcp,
                           enloc_sqtl_max_rcp_intron = combined_intron_names,
                           enloc_sqtl_max_rcp_lead_coloc_SNP = combined_lead_coloc_SNP_names,
                           enloc_sqtl_max_rcp_lead_snp_rcp = combined_lead_snp_rcp_names)
    } else if (nrow(best_df) == 0){
      best_final <- tibble(enloc_sqtl_max_rcp_tissue = NA,
                           enloc_sqtl_max_locus_rcp = NA,
                           enloc_sqtl_min_locus_rcp = NA,
                           enloc_sqtl_max_rcp_intron = NA,
                           enloc_sqtl_max_rcp_lead_coloc_SNP = NA,
                           enloc_sqtl_max_rcp_lead_snp_rcp = NA)
    } else {
      best_df$lead_snp_rcp <- as.character(best_df$lead_snp_rcp)
      best_df$enloc_sqtl_min_locus_rcp = min_rcp
      best_final <- best_df %>% rename(enloc_sqtl_max_rcp_tissue = tissue, 
                                       enloc_sqtl_max_locus_rcp = locus_rcp,
                                       enloc_sqtl_max_rcp_lead_coloc_SNP = lead_coloc_SNP,
                                       enloc_sqtl_max_rcp_lead_snp_rcp = lead_snp_rcp)
    }
    final <- bind_cols(breast_final, best_final)
    # browser()
    return(final)
  }

  # debug_df <- c("ENSG00000164440.14", "ENSG00000226571.1", "ENSG00000265531.3", "ENSG00000166924.8") %>%
  #   map_dfr(~get_max_tissue_sqtl(., enloc_sqtl_df, tissue_filter=tissue_filter))
  # browser()
  t_out_df <- unique(enloc_sqtl_df$group) %>% 
    map_dfr(~get_max_tissue_sqtl(., enloc_sqtl_df, tissue_filter=tissue_filter))

  opath <- file.path(odir, glue("gene_rolled_enloc_sqtl_{fname_only}"))
  write_delim(t_out_df, opath, delim="\t")
  return(t_out_df)
}

if (!interactive()) {
  library(purrr)
  ENLOC_RESULTS %>% walk(~main(.))
} else {
  gimme <- main(ENLOC_RESULTS[4])
  # og <- fread(ENLOC_RESULTS[1]) %>%
  #   rename(intron = molecular_qtl_trait) %>%
  #   as.tibble()
}
