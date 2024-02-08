ENLOC_RESULTS <- list.files(path="../../03_enloc/output/final_enloc_result_eqtl", pattern="*.txt",
                            full.names=TRUE)

library(tidyverse)
library(data.table)
library(glue)

main <- function(input_enloc_eqtl_path,
                 tissue_filter="../../02_acat_eqtl_sqtl/input/tissue_lists/11tiss.txt",
                 odir="../output"){
  require(tidyverse)
  require(data.table)
  require(glue)
  require(purrr)

  print(input_enloc_eqtl_path)
  fname_only <- basename(input_enloc_eqtl_path)

  enloc_eqtl_df <- fread(input_enloc_eqtl_path) %>%
    rename(group = molecular_qtl_trait) %>%
    as.tibble()

  get_max_tissue_eqtl <- function(gname, enloc_eqtl, enloc_eqtl_breast, tissue_filter=tissue_filter){
    # if (gname == "ENSG00000198807.12"){
    #   browser()
    # }
    if (!is.na(tissue_filter)){
      tiss_filter <- fread(tissue_filter, header=F)
      enloc_eqtl <- enloc_eqtl %>% 
        filter(tissue %in% tiss_filter$V1)
    }
    rdf <- enloc_eqtl %>% filter(group == gname) %>%
      filter(locus_rcp == max(locus_rcp, na.rm=T)) %>%
      select(group, tissue, locus_rcp, lead_coloc_SNP, lead_snp_rcp)

    min_rcp <- rdf %>% filter(group == gname) %>%
      pull(locus_rcp) %>% min(., na.rm=T)
  
    if (nrow(rdf) > 1){
      combined_tiss_names <- paste(unique(rdf$tissue) %>% sort(), collapse=", ")
      combined_lead_coloc_SNP_names <- paste(unique(rdf$lead_coloc_SNP) %>% sort(), collapse=", ")
      combined_lead_snp_rcp_names <- paste(unique(rdf$lead_snp_rcp) %>% sort(), collapse=", ")
      final <- tibble(group = gname, enloc_eqtl_max_rcp_tissue = combined_tiss_names,
                      enloc_eqtl_lead_coloc_SNP = combined_lead_coloc_SNP_names,
                      enloc_eqtl_lead_snp_rcp = combined_lead_snp_rcp_names,
                      enloc_eqtl_max_locus_rcp = rdf$locus_rcp[1],
                      enloc_eqtl_min_locus_rcp = min_rcp)
    } else if (nrow(rdf) == 0){
      final <- tibble(group = gname, enloc_eqtl_max_rcp_tissue = NA,
                      enloc_eqtl_lead_coloc_SNP = NA,
                      enloc_eqtl_lead_snp_rcp = NA,
                      enloc_eqtl_max_locus_rcp = NA,
                      enloc_eqtl_min_locus_rcp = NA)
    } else {
      rdf$lead_snp_rcp <- as.character(rdf$lead_snp_rcp)
      rdf$enloc_eqtl_min_locus_rcp = min_rcp
      final <- rdf %>% rename(enloc_eqtl_max_rcp_tissue = tissue, enloc_eqtl_max_locus_rcp = locus_rcp,
                              enloc_eqtl_lead_coloc_SNP = lead_coloc_SNP,
                              enloc_eqtl_lead_snp_rcp = lead_snp_rcp)
    }

    breast_portion <- enloc_eqtl_breast %>%
      filter(group == gname) %>%
      rename(enloc_eqtl_breast_locus_rcp = locus_rcp,
             enloc_eqtl_breast_lead_coloc_SNP = lead_coloc_SNP,
             enloc_eqtl_breast_lead_snp_rcp = lead_snp_rcp) %>%
      select(-c("group"))

    breast_portion <- breast_portion %>%
      filter(enloc_eqtl_breast_locus_rcp == max(enloc_eqtl_breast_locus_rcp, na.rm=T))

    if (nrow(breast_portion) > 1){
      combined_lead_coloc_SNP_names <- paste(unique(breast_portion$enloc_eqtl_breast_lead_coloc_SNP) %>% sort(), collapse=", ")
      combined_lead_snp_rcp_names <- paste(unique(breast_portion$enloc_eqtl_breast_lead_snp_rcp) %>% sort(), collapse=", ")
      breast_final <- tibble(
                      enloc_eqtl_breast_lead_coloc_SNP = combined_lead_coloc_SNP_names,
                      enloc_eqtl_breast_lead_snp_rcp = combined_lead_snp_rcp_names,
                      enloc_eqtl_breast_max_locus_rcp = breast_portion$enloc_eqtl_breast_locus_rcp[1])
    } else if(nrow(breast_portion) == 0){
      breast_final <- tibble(
                      enloc_eqtl_breast_lead_coloc_SNP = NA,
                      enloc_eqtl_breast_lead_snp_rcp = NA,
                      enloc_eqtl_breast_max_locus_rcp = NA)
    } else {
      breast_portion$enloc_eqtl_breast_lead_snp_rcp <- as.character(breast_portion$enloc_eqtl_breast_lead_snp_rcp)
      breast_final <- breast_portion %>% rename(enloc_eqtl_breast_max_locus_rcp = enloc_eqtl_breast_locus_rcp,
                              enloc_eqtl_breast_lead_coloc_SNP = enloc_eqtl_breast_lead_coloc_SNP,
                              enloc_eqtl_breast_lead_snp_rcp = enloc_eqtl_breast_lead_snp_rcp)
    }

    return(bind_cols(final, breast_final))
  }

  eqtl_breast_df <- enloc_eqtl_df %>%
    filter(tissue == "Breast_Mammary_Tissue")
  eqtl_best_locus_rcp_tissue <- unique(enloc_eqtl_df$group) %>% 
    map_dfr(~get_max_tissue_eqtl(., enloc_eqtl_df, enloc_eqtl_breast=eqtl_breast_df, tissue_filter=tissue_filter))

  t_out_df <- eqtl_best_locus_rcp_tissue

  opath <- file.path(odir, glue("gene_rolled_enloc_eqtl_{fname_only}"))
  write_delim(t_out_df, opath, delim="\t")
  return(list(final=eqtl_best_locus_rcp_tissue, breast=eqtl_breast_df))
}

if (!interactive()) {
  library(purrr)
  ENLOC_RESULTS %>% walk(~main(.))
} else {
  gimme <- main(ENLOC_RESULTS[1])
  # og <- fread(ENLOC_RESULTS[1]) %>%
  #   rename(group = molecular_qtl_trait) %>%
  #   as.tibble()
}
