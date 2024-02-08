main <- function(known_vars="../input/GWAS_variant_annotated_James_DH_2023Jan5.tsv",
                 ofile="../output/gencode_26_40_conversion_and_distances_to_known_brca_index_snps.tsv"){
  require(tidyverse)
  require(data.table)
  require(glue)

  gc_26_df <- fread("../input/gencode_v26_all.txt") %>%
    mutate(gene_id_no_decimal = str_extract(gene_id, "(ENSG\\d{1,})"))
  # chromosome,source,feature_type,start_location,end_location,gene_id,gene_type,gene_name
  gc_38_df <- fread("../input/gencode_v38_all.txt") %>%
    mutate(gene_id_no_decimal = str_extract(gene_id, "(ENSG\\d{1,})")) %>%
    rename(chromosome = seqname)
  gc_40_df <- fread("../input/gencode_v40_all.txt") %>%
    mutate(gene_id_no_decimal = str_extract(gene_id, "(ENSG\\d{1,})"))

  gc_combined_df <- left_join(gc_26_df, gc_40_df, by = c("gene_id_no_decimal", "chromosome"),
                              suffix=c(".v26", ".v40")) %>%
    mutate(gencode_conversion_success = ifelse(is.na(gene_id.v40), FALSE, TRUE))
  gc_combined_38_df <- left_join(gc_26_df %>% select("gene_id_no_decimal", "chromosome", "gene_name"),
                                 gc_38_df %>% select("gene_id_no_decimal", "chromosome", "gene_name"),
                                 by = c("gene_id_no_decimal", "chromosome"),
                              suffix=c(".v26", ".v38"))
  gc_combined_df <- left_join(gc_combined_df, gc_combined_38_df %>% select("gene_id_no_decimal", "chromosome", "gene_name.v38"),
                              by = c("gene_id_no_decimal", "chromosome"))

  # Chromosome 'hg38 position'
  brca_known_vars_df <- fread(known_vars, fill=T)


  # Add Cytoband locus
  e <- data.table()
  cytoband <- fread("../input/cytoBand.txt.gz", col.names=c("Chr", "start_position", "end_location", "locus", "info"))
  for (i in 1:22){
    a <- gc_combined_df %>% filter(chromosome == glue("chr{i}")) %>%
      mutate(start_location = ifelse(is.na(`start_location.v40`), `start_location.v26`, `start_location.v40`),
             end_location = ifelse(is.na(`end_location.v40`), `end_location.v26`, `end_location.v40`))
    b <- cytoband %>% filter(Chr == glue("chr{i}"))
    for (j in 1:nrow(b)){
      start <- b[j,]$start_position
      stop <- b[j,]$end_location
      d <- a %>% filter(start_location>= start & start_location <= stop) %>%
        mutate(`Gene Loci` = paste(i,b[j,]$locus, sep=""))
      e <- rbind(e, d)
    }
  }
  gc_combined_df <- e

  # Get closest index SNP
  f <- data.table()
  g <- data.table()
  # for (chr_num in 1:1){
  for (chr_num in 1:22){
    a <- gc_combined_df %>% filter(chromosome == glue("chr{chr_num}")) %>%
      mutate(start_location = ifelse(is.na(`start_location.v40`), `start_location.v26`, `start_location.v40`),
             end_location = ifelse(is.na(`end_location.v40`), `end_location.v26`, `end_location.v40`),
             gene_name = ifelse(is.na(`gene_name.v40`), `gene_name.v26`, `gene_name.v40`)
             )
    b <- brca_known_vars_df %>% filter(Chromosome == chr_num)

    # for (j in 1:1){
    for (j in 1:nrow(a)){
      f$t_start_s <- a$start_location[j]
      f$t_stop_s <- a$end_location[j]

      f$gene <- a$gene_id_no_decimal[j] 

      c <- b %>% mutate(distance_from_start = f$t_start_s[1] - `hg38 position`,
                        distance_from_stop = f$t_stop_s[1] - `hg38 position`)
      c <- c %>% mutate(abs_distance_from_start = abs(distance_from_start),
                        abs_distance_from_stop = abs(distance_from_stop))
      c <- c %>% filter(abs_distance_from_start == min(abs_distance_from_start))

      c2 <- b %>% mutate(distance_from_stop = f$t_stop_s[1] - `hg38 position`,
                         distance_from_start = f$t_start_s[1] - `hg38 position`)
      c2 <- c2 %>% mutate(abs_distance_from_stop = abs(distance_from_stop),
                          abs_distance_from_start = abs(distance_from_start))
      c2 <- c2 %>% filter(abs_distance_from_stop == min(abs_distance_from_stop))

      if (nrow(c) > 1){
        # print(c)
        c <- c %>% filter(`First Author` == "Jia G")
      }
      if (nrow(c2) > 1){
        # print(c2)
        c2 <- c2 %>% filter(`First Author` == "Jia G")
      }

      if (c$abs_distance_from_start <= c2$abs_distance_from_stop){
        closest_snp <- c
      } else {
        closest_snp <- c2
      }

      f$SNP <- closest_snp$`Index SNP`
      f$`Index SNP Loci` <- closest_snp$Loci

      # Case 1, SNP inside gene:
      #   (distance from start negative, from stop, postive) distance_to_snp = 0
      # Case 2, SNP to right of gene:
      #   (distance from start negative, from stop, negative) distance_to_snp = abs(distance_from_stop)
      # Case 3, SNP to left of gene:
      #   (distance from start negative, from stop, postive) distance_to_snp = abs(distance_from_start)
      if ((closest_snp$distance_from_start <= 0) & (closest_snp$distance_from_stop >= 0)){
        f$distance_to_snp <- 0
      } else if ((closest_snp$distance_from_start <= 0) & (closest_snp$distance_from_stop <= 0)){
        f$distance_to_snp <- closest_snp$abs_distance_from_stop
      } else if ((closest_snp$distance_from_start >= 0) & (closest_snp$distance_from_stop >= 0)){
        f$distance_to_snp <- closest_snp$abs_distance_from_start
      } else {
        print(closest_snp)
        print("Investigate further.")
        browser()
      }
      f$distance_to_start <- closest_snp$abs_distance_from_start
      f$distance_to_stop <- closest_snp$abs_distance_from_stop
      f <- data.table(f$t_start_s, f$t_stop_s, f$gene, f$SNP, f$`Index SNP Loci`, f$distance_to_snp, f$distance_to_start, f$distance_to_stop)
      g <- bind_rows(g,f)
    }
  }
  colnames(g)<- c("t_start_s","t_stop_s","gene","Index SNP", "Index SNP Loci", "distance_to_snp", "distance_to_start", "distance_to_stop")

  # gc_combined_df <- left_join(gc_combined_df, g, by = c("gene_id_no_decimal" = "gene")) %>%
  # Inner join to leave out chromosomes that aren't 1-22
  gc_combined_df <- inner_join(gc_combined_df, g, by = c("gene_id_no_decimal" = "gene")) %>%
    select(-t_start_s, -t_stop_s) %>%
    mutate(`Distance from Index SNP(kb)` = distance_to_snp / 1000,
           `Distance from Gene Start(kb)` = distance_to_start / 1000,
           `Distance from Gene End(kb)` = distance_to_stop / 1000)
  prev_twas_df <- fread("../input/Ten_TWAS_papers_gene_list_2023Feb28.tsv")

  gc_combined_df <- gc_combined_df %>%
    mutate(brca_twas_reported = ifelse(gene_id_no_decimal %in% prev_twas_df$ENSG, T, F))
  write_delim(gc_combined_df, ofile, delim="\t")
  return(gc_combined_df)
}


if (!interactive()) {
  main("../input/bcac_ukb_meta_add_new_vars_GWAS_variant_annotated_James_DH_2023Jun02.tsv",
       ofile="../output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps_2023Jun02.tsv")
  # main("../input/GWAS_variant_annotated_James_DH_2023Jun02.tsv")
} else {
  main("../input/bcac_ukb_meta_add_new_vars_GWAS_variant_annotated_James_DH_2023Jun02.tsv",
       ofile="../output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps_2023Jun02.tsv")
  # main("../input/GWAS_variant_annotated_James_DH_2023Jun02.tsv")
}
