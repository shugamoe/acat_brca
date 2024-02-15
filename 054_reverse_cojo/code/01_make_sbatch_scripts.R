if (!exists("ERPOS_MA_DF")){
	ERPOS_MA_DF <- data.table::fread("../input/BCAC_ERPOS_BreastCancer_EUR_gwas.ma")
}
if (!exists("ERNEG_MA_DF")){
	ERNEG_MA_DF <- data.table::fread("../input/meta_analysis_BCAC_CIMBA_erneg_brca_gwas.ma")
}

main <- function(tiss_path="../input/tissue_lists/11tiss.txt"){
	require(tidyverse)
	require(data.table)
	require(readxl)
	require(glue)

	tiss_list_basename_no_ext <- basename(tiss_path) %>%
	  str_replace_all("\\.txt", "")
	tiss_vec <- fread(tiss_path, header=F)$V1

	if (!exists("ERPOS_MA_DF")){
	  erpos_ma_df <- fread("../input/BCAC_ERPOS_BreastCancer_EUR_gwas.ma")
	} else {
	  erpos_ma_df <- ERPOS_MA_DF
	}

	if (!exists("ERNEG_MA_DF")){
	  erneg_ma_df <- fread("../input/meta_analysis_BCAC_CIMBA_erneg_brca_gwas.ma")
	} else {
	  erneg_ma_df <- ERNEG_MA_DF
	}

	rsid_convert_df <- fread("../input/og_cojo_input/bcac_ukb_meta_add_new_vars_GWAS_variant_annotated_James_DH_2023Jun02.tsv",
													 fill=T) %>%
	  select(`Chromosome`, `hg38 position`, `Index SNP`)

	erpos_master <- read_excel("../input/Additional_File_1_2024Feb13.xlsx", sheet="TABLE S3", range="A2:G791") %>%
		filter(`Region of TWAS or GWAS signals` == "Both") %>%
		mutate(chr_num = str_extract(Locus, "^\\d{1,2}")) %>%
		mutate(chromosome = glue("chr{chr_num}")) %>%
		distinct() %>%
		left_join(rsid_convert_df, by = c("GWAS SNP" = "Index SNP")) %>%
		mutate(chromosome_position=glue("chr{Chromosome}_{`hg38 position`}"),
					 Locus = str_replace(Locus, " ", "_"))

	erneg_master <- read_excel("../input/Additional_File_1_2024Feb13.xlsx", sheet="TABLE S4", range="A2:G251") %>%
		filter(`Region of TWAS or GWAS signals` == "Both") %>%
		mutate(chr_num = str_extract(Locus, "^\\d{1,2}")) %>%
		mutate(chromosome = glue("chr{chr_num}")) %>%
		distinct() %>%
		left_join(rsid_convert_df, by = c("GWAS SNP" = "Index SNP")) %>%
		mutate(chromosome_position=glue("chr{Chromosome}_{`hg38 position`}"),
					 Locus = str_replace(Locus, " ", "_"))

	write_tsv(erpos_master, "../output/ERPOS_master.tsv")
	write_tsv(erneg_master, "../output/ERNEG_master.tsv")

	# Get the ensg v26 ids for the genes
	# Lots of repeated code from 0000_get_ensg_v26_ids_from_S3andS4.R, but in a rush
	og_erpos_genes <- erpos_master %>% filter(!is.na(`TWAS gene`)) %>% nrow()
	og_erneg_genes <- erneg_master %>% filter(!is.na(`TWAS gene`)) %>% nrow()

	gname_ensg_df <- fread("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/999_reporting_friendly_tables/output/bcac_ukb_meta_add_new_vars_gencode_26_40_conversion_and_distances_to_known_brca_index_snps.tsv") %>%
		select(`chromosome`, `gene_id.v26`, `gene_name.v40`, `gene_name.v38`, `gene_name.v26`) %>%
		distinct()

	erpos_master_conv <- erpos_master %>%
		filter(!is.na(`TWAS gene`)) %>%
		left_join(gname_ensg_df %>% select(`chromosome`, `gene_name.v40`, `gene_id.v26`), by = c("TWAS gene" = "gene_name.v40", "chromosome")) %>% 
		filter(!is.na(`gene_id.v26`))

	erpos_master_conv2 <- erpos_master %>%
		filter(!is.na(`TWAS gene`)) %>%
		left_join(gname_ensg_df %>% select(`chromosome`, `gene_name.v26`, `gene_id.v26`), by = c("TWAS gene" = "gene_name.v26", "chromosome")) %>% 
		filter(!is.na(`gene_id.v26`))

	erpos_master_all_conv <- bind_rows(erpos_master_conv, erpos_master_conv2) %>%
		distinct()

	erpos_all_conv_num <- nrow(erpos_master_all_conv)

	if (og_erpos_genes != erpos_all_conv_num){
		browser()
		error("Some erpos genes not converted")
	}

	erpos_master <- erpos_master %>%
		left_join(erpos_master_all_conv %>% select(`TWAS gene`, `gene_id.v26`))

	erneg_master_conv <- erneg_master %>%
		filter(!is.na(`TWAS gene`)) %>%
		left_join(gname_ensg_df %>% select(`chromosome`, `gene_name.v40`, `gene_id.v26`), by = c("TWAS gene" = "gene_name.v40", "chromosome")) %>% 
		filter(!is.na(`gene_id.v26`))

	erneg_master_conv2 <- erneg_master %>%
		filter(!is.na(`TWAS gene`)) %>%
		left_join(gname_ensg_df %>% select(`chromosome`, `gene_name.v26`, `gene_id.v26`), by = c("TWAS gene" = "gene_name.v26", "chromosome")) %>% 
		filter(!is.na(`gene_id.v26`))

	erneg_master_all_conv <- bind_rows(erneg_master_conv, erneg_master_conv2) %>%
		distinct()

	erneg_all_conv_num <- nrow(erneg_master_all_conv)

	if (og_erneg_genes != erneg_all_conv_num){
		browser()
		error("Some erneg genes not converted")
	}

	erneg_master <- erneg_master %>%
		left_join(erneg_master_all_conv %>% select(`TWAS gene`, `gene_id.v26`))

	process_master <- function(cur_master, cur_ma_df, gwas_type){
		uniq_loci <- cur_master$`Locus` %>% unique()

		if (gwas_type == "ERPOS"){
			ma_path <- "../input/BCAC_ERPOS_BreastCancer_EUR_gwas.ma"
		} else if (gwas_type == "ERNEG"){
			ma_path <- "../input/meta_analysis_BCAC_CIMBA_erneg_brca_gwas.ma"
		}

		for (cur_loci in uniq_loci){
			cur_loci_master <- cur_master %>%
				filter(Locus == cur_loci)
			cur_loci_chrom <- cur_loci_master %>%
				filter(!is.na(`Chromosome`)) %>%
				head(1) %>%
				pull(`Chromosome`) 

			cur_loci_twas_genes <- cur_loci_master %>%
				filter(!is.na(`TWAS gene`)) %>%
				pull(`gene_id.v26`) %>% 
				unique()
			for (model_type in c("eqtl", "sqtl", "both")){ # Run COJO with the predictive SNPs from eqtl only, sqtl only, or both model types
				index_snps <- get_index_snps(cur_loci_twas_genes, tiss_vec, cur_ma_df, gwas_type, model_type)  
				if (nrow(index_snps) == 0){
					message(glue("No valid prediction SNPs for genes in {gwas_type} | {cur_loci} | {model_type} |{tiss_list_basename_no_ext}"))
					next
				}
				want_snps <- cur_loci_master %>% # Want conditional effects to be calced for these SNPs
					filter(!is.na(`GWAS SNP`)) %>%
					select(`chromosome_position`)

				extract_snps <- bind_rows(index_snps, want_snps) %>%
					select(`chromosome_position`) %>%
					distinct()

				# Write some files for later COJO/clumping/COJO with clumped data
				gp_id <- glue("{tiss_list_basename_no_ext}__{gwas_type}__{cur_loci}__gtex_{model_type}")
				want_snp_path <- glue("../output/intermediate_data/raw_cojo/{gp_id}.want") 
				extract_snp_path <- glue("../output/intermediate_data/raw_cojo/{gp_id}.extract") 
				index_snp_clump_path <- glue("../output/intermediate_data/plink_clumping/{gp_id}.snps.clump.target")
				index_snp_cojo_path <- glue("../output/intermediate_data/raw_cojo/{gp_id}.snps.index")
				write_tsv(want_snps, want_snp_path,
									col_names=F)
				write_tsv(extract_snps, extract_snp_path,
									col_names=F)
				write_delim(index_snps, index_snp_clump_path, delim=" ")
				write_tsv(index_snps %>% select(`chromosome_position`), index_snp_cojo_path)

				olap_df <- inner_join(want_snps, index_snps)
				if (nrow(olap_df) > 0){
					message(glue("Index/Cond SNP overlap in {gwas_type} | {cur_loci} | {model_type} | {tiss_list_basename_no_ext}"))
					write_delim(olap_df %>% select(`chromosome_position`), glue("../output/intermediate_data/overlap_tracking/{gp_id}.tsv"))
				}

				# With want_snps and index snps in hand, write a file for running COJO with no clumping conditional list
				raw_cojo_command <- glue("
				#!/bin/bash
				#SBATCH --job-name={gp_id}
				#SBATCH --time=0:30:00
				#SBATCH --nodes=1
				#SBATCH --cpus-per-task=1
				#SBATCH --tasks-per-node=1
				#SBATCH --mem=16GB
				#SBATCH --output=logs/raw_cojo/%x.out

				cd $SLURM_SUBMIT_DIR

				module load gcc/12.1.0
				module load miniconda3/23.1.0
				source activate /gpfs/data/huo-lab/jmcclellan/software/envs/parsl
				export PATH=/gpfs/data/gao-lab/Julian/software/gcta-1.94.1-linux-kernel-3-x86_64:$PATH

				gcta-1.94.1 \\
				  --bfile ../input/og_cojo_input/backup_COJO_REF_LD/hg38_EUR_chr_{cur_loci_chrom}_gtex_v8_1000G.phase3.genotypes.final \\
				  --cojo-file {ma_path} \\
					--cojo-cond {index_snp_cojo_path} \\
					--extract {extract_snp_path} \\
					--out ../output/intermediate_data/raw_cojo/{gp_id}
				")
				writeLines(raw_cojo_command, glue("jobs/raw_cojo/{gp_id}.sbatch"))

				for (want_r2 in c(.65, .5, .4, .3, .2, .1)){
				  gp_id <- glue("{tiss_list_basename_no_ext}__{gwas_type}__{cur_loci}__gtex_{model_type}__clump_r2_{want_r2}")
					# Write a file that will do plink clumping on the index SNPs
					clump_prefix <- glue("../output/intermediate_data/plink_clumping/{gp_id}")
					expected_clump_output <- glue("{clump_prefix}.clumped")
					plink_clump_command <- glue("
					#!/bin/bash
					#SBATCH --job-name={gp_id}
					#SBATCH --time=0:30:00
					#SBATCH --nodes=1
					#SBATCH --cpus-per-task=1
					#SBATCH --tasks-per-node=1
					#SBATCH --mem=16GB
					#SBATCH --output=logs/plink_clumping/%x.out

					cd $SLURM_SUBMIT_DIR

					module load gcc/12.1.0
					module load miniconda3/23.1.0
					module load plink/1.9


					# Clump using same settings as with original COJO, BRCA
					plink --bfile ../input/og_cojo_input/backup_COJO_REF_LD/hg38_EUR_chr_{cur_loci_chrom}_gtex_v8_1000G.phase3.genotypes.final \\
						--clump {index_snp_clump_path} \\
						--clump-r2 {want_r2} \\
						--clump-kb 10001 \\
						--clump-p1 1 \\
						--clump-p2 1 \\
						--memory 16000 \\
						--threads 1 \\
						--clump-snp-field chromosome_position \\
						--clump-field P \\
						--out {clump_prefix}

					if test -f {expected_clump_output}; then
						tail -n +2 {expected_clump_output} | tr -s ' ' | cut -d ' ' -f 4 | grep 'chr' > {expected_clump_output}.trimmed
					else
						echo 'No plink output?' # Use the most significant old conditioned SNP instead I guess
						cat {index_snp_clump_path} | tail -n +2 | cut -d ' ' -f 1 | tail -1  > {expected_clump_output}.trimmed
					fi
					")
					writeLines(plink_clump_command, glue("jobs/plink_clumping/{gp_id}.sbatch"))

					# Write a file that will run COJO based off of clumped index SNPs
					clump_extract_path <- glue("../output/intermediate_data/clumped_cojo/{gp_id}.extract")
					clumped_cojo_command <- glue("
					#!/bin/bash
					#SBATCH --job-name={gp_id}
					#SBATCH --time=0:30:00
					#SBATCH --nodes=1
					#SBATCH --cpus-per-task=1
					#SBATCH --tasks-per-node=1
					#SBATCH --mem=16GB
					#SBATCH --output=logs/clumped_cojo/%x.out

					cd $SLURM_SUBMIT_DIR

					cat {expected_clump_output}.trimmed {want_snp_path} | sort --unique > {clump_extract_path}

					module load gcc/12.1.0
					module load miniconda3/23.1.0
					source activate /gpfs/data/huo-lab/jmcclellan/software/envs/parsl
					export PATH=/gpfs/data/gao-lab/Julian/software/gcta-1.94.1-linux-kernel-3-x86_64:$PATH

					gcta-1.94.1 \\
						--bfile ../input/og_cojo_input/backup_COJO_REF_LD/hg38_EUR_chr_{cur_loci_chrom}_gtex_v8_1000G.phase3.genotypes.final \\
						--cojo-file {ma_path} \\
						--cojo-cond {expected_clump_output}.trimmed \\
						--extract {clump_extract_path} \\
						--out ../output/intermediate_data/clumped_cojo/{gp_id}
					")
					writeLines(clumped_cojo_command, glue("jobs/clumped_cojo/{gp_id}.sbatch"))
				}
			}
		}
	}

	process_master(erpos_master, erpos_ma_df, "ERPOS")
	process_master(erneg_master, erneg_ma_df, "ERNEG")
}


get_index_snps <- function(ensg_v26_ids, tiss_vec, cojo_ma_df, gwas_type, model_type){
	return_df <- tibble()

  if (model_type == "both"){
		model_types <- c("eqtl", "sqtl")
	} else {
		model_types <- model_type
	}

	for (mtype in model_types){
		for (ctiss in tiss_vec){
			for (cgene in ensg_v26_ids){
				check_path <- glue("../output/intermediate_data/gene_snp_lists/{gwas_type}/{mtype}/{ctiss}/{cgene}.snplist.all.for.prediction")
				if (file.exists(check_path)){
					cur_df <- fread(check_path, header=F)
					return_df <- bind_rows(return_df, cur_df)
				}
			}
		}
	}

	if (!("V1" %in% names(return_df))){
		return(tibble())
	}
	return_df <- return_df %>%
		distinct() %>%
		rename(`chromosome_position` = `V1`)

	final_return_df <- return_df %>%
		inner_join(cojo_ma_df) %>%
		select(`chromosome_position`, `pvalue`) %>%
		rename(`P` = `pvalue`)

	return(final_return_df)
}

if (interactive()){
	main()
} else {
	main()
}
