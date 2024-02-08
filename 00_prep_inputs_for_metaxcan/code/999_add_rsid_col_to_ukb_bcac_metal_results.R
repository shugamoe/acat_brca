YOUR_METAL_FILE="/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/00_prep_inputs_for_metaxcan/output/metal/final_metal_bcac-white_ukb1.txt"
OUT_METAL_PATH="/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/00_prep_inputs_for_metaxcan/output/metal/rsid_final_metal_bcac-white_ukb1.txt"

main <- function(metal_file_path, out_file_path,
                 metal_link_col="MarkerName",
                 rsid_source_link_col="var_name"
                 ){
  require(tidyverse)
  require(data.table)
  require(glue)
  
  if (metal_file_path==""){
    print("Set YOUR_METAL_FILE at the top of the file, or run in an R interpreter and set <metal_file_path>")
    return()
  }
  if (out_file_path==""){
    print("Set OUT_METAL_FILE at the top of the file, or run in an R interpreter and set <out_file_path>")
    return()
  }

  # MarkerName  Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Effect  StdErr  P-value Direction   HetISq  HetChiSq    HetDf   HetPVal
  metal_df <- fread(metal_file_path)
  print(glue("Read file: {metal_file_path}"))

  # CHROM POS rsID altID linkID REF ALT MinorAllele MAF Info A1 OBS_CT OR beta_ALT SE Z_STAT P var_name EAF 
  ukb_path <- "/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/00_prep_inputs_for_metaxcan/input/Anno_res_anyBC_full_chrom.breast.logistic.gz"

  ukb_rsid_key <- fread(ukb_path, select=c("rsID", "var_name")) %>%
    filter(str_detect(rsID, "rs"))
  print(glue("Read file: {ukb_path}"))

  bcac_path <- "/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/00_prep_inputs_for_metaxcan/input/oncoarray_bcac_public_release_oct17.rsid.conv.txt"
  bcac_rsid_key <- fread(bcac_path,
                         select=c("var_name", "phase3_1kg_id")) %>%
    filter(str_detect(phase3_1kg_id, "rs")) %>%
    separate(phase3_1kg_id, c("rsID", NA, NA, NA), sep=":", remove=F) %>%
    select(all_of(c("rsID", "var_name")))
  print(glue("Read file: {bcac_path}"))

  master_rsid_key <- bind_rows(ukb_rsid_key, bcac_rsid_key) %>%
    distinct()

                 # metal_link_col="MarkerName",
                 # rsid_source_link_col="var_name"
  metal_df_with_rsid = left_join(metal_df, master_rsid_key,
                                 # https://stackoverflow.com/questions/28399065/dplyr-join-on-by-a-b-where-a-and-b-are-variables-containing-strings
                                 by = setNames(nm=metal_link_col, rsid_source_link_col))

  final_col_order <- c("MarkerName", "rsID", "Allele1", "Allele2", "Freq1", 
                       "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr",
                       "P-value", "Direction",   "HetISq",  "HetChiSq",
                       "HetDf", "HetPVal")
  write_delim(metal_df_with_rsid %>% select(all_of(final_col_order)), out_file_path, delim="\t")
  print(glue("File with RSID written at {out_file_path}"))
}

if (!interactive()) {
  main(metal_file_path=YOUR_METAL_FILE,
       out_file_path=OUT_METAL_PATH
       )
}
