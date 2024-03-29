python3 -m pdb /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/conditional_covariance_for_model_group.py \
        -parquet_genotype_folder ../../00_prep_inputs_for_metaxcan/input/gtex_v8_parquet_eur_maf0.01_biallelic/ \
        -parquet_genotype_pattern "gtex_v8_eur_itm.chr(\d+).variants.parquet" \
        -group ../input/combine_phenotype_groups_multi_tiss.txt.gz \
        -model_db_group_key ../input/combine_eqtl_for_covar_calc.db \
        -model_db_group_value ../input/combine_sqtl_for_covar_calc.db \
        -parsimony 8 \
        -want_genes ../output/intermediate_data/condTWAS_cond_covar/meta_analysis_BCAC_UKB_ovr_brca/want_genes.txt \
        -condition_info_dir ../output/intermediate_data/condTWAS_cond_covar/meta_analysis_BCAC_UKB_ovr_brca \
        -pred_pattern {want_gene}.snps.pred.txt \
        -cond_pattern {want_gene}.snps.condition.txt \
        --individuals ../../02_acat_eqtl_sqtl/input/gtex_v8_by_tiss_individuals/Adipose_Visceral_Omentum_individuals.txt \
        -tissue Adipose_Visceral_Omentum \
        -output_dir tmp/
