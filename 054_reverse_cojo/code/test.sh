	# --model_db_path ../input/gtex_v8_sqtl_dbs_mashr/mashr_${tiss}.db \
	# --model_db_path ../input/gtex_v8_sqtl_dbs_mashr/mashr_${tiss}.db \
tiss="Whole_Blood"
	python3 /gpfs/data/huo-lab/jmcclellan/software/MetaXcan/software/SMulTi_gene_snps.py \
	    --gwas_file /gpfs/data/huo-lab/jmcclellan/acat_brca_combined/00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz \
		--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
		--grouping /gpfs/data/gao-lab/Julian/gaolab_hub/projects/gtex_covar_extension/alvaro/single_tiss_prep_for_multi_tiss_covar_calc/${tiss}.leafcutter.phenotype_groups.txt.gz GTEx_sQTL \
		--model_db_path ../input/gtex_v8_sqtl_dbs_mashr/mashr_${tiss}.db \
		--keep_non_rsid --model_db_snp_key varID \
		--associations /gpfs/data/huo-lab/jmcclellan/acat_brca_combined/01_spredixcan_eqtl_sqtl/output/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_meta_analysis_BCAC_CIMBA_erneg_brca__PM__${tiss}.csv "SPrediXcan" \
		--verbosity 10 \
		--snp_list_prefix ../output/intermediate_data/gene_snp_lists/ERNEG/sqtl/${tiss} \
		--output ../output/intermediate_data/gene_snp_lists/ERNEG/sqtl/${tiss}/${tiss} \
		--gene_whitelist ../input/ERNEG_ensg_v26_ids.txt "gene_id.v26"
