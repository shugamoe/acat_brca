#!/bin/bash
#SBATCH --job-name=ERNEG_get_snplists
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=16GB
#SBATCH --output=logs/%x.out

cd $SLURM_SUBMIT_DIR

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/huo-lab/jmcclellan/software/envs/parsl
	    # --model_db /gpfs/data/huo-lab/jmcclellan/acat_brca_combined/04_cojo_ctwas/input/combine_eqtl_for_covar_calc.db \

for tiss in `cat ../input/tissue_lists/all_tissues.txt`; do
	echo "get SQTL"
	python3 /gpfs/data/huo-lab/jmcclellan/software/MetaXcan/software/SMulTi_gene_snps.py \
	    --gwas_file /gpfs/data/huo-lab/jmcclellan/acat_brca_combined/00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz \
		--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
		--grouping /gpfs/data/gao-lab/Julian/gaolab_hub/projects/gtex_covar_extension/alvaro/single_tiss_prep_for_multi_tiss_covar_calc/${tiss}.leafcutter.phenotype_groups.txt.gz GTEx_sQTL \
		--model_db_path ../input/gtex_v8_sqtl_dbs_mashr/mashr_${tiss}.db \
		--keep_non_rsid --model_db_snp_key varID \
		--associations /gpfs/data/huo-lab/jmcclellan/acat_brca_combined/01_spredixcan_eqtl_sqtl/output/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_meta_analysis_BCAC_CIMBA_erneg_brca__PM__${tiss}.csv "SPrediXcan" \
		--verbosity 10 \
		--snp_list_prefix "../output/intermediate_data/gene_snp_lists/ERNEG/sqtl/${tiss}/" \
		--output ../output/intermediate_data/gene_snp_lists/ERNEG/sqtl/${tiss}/${tiss} \
		--gene_whitelist ../input/ERNEG_ensg_v26_ids.txt "gene_id.v26"

	echo "get EQTL"
	python3 /gpfs/data/huo-lab/jmcclellan/software/summary-gwas-imputation/src/groups_and_conditioned_covariance_for_model.py \
	    --gwas_file /gpfs/data/huo-lab/jmcclellan/acat_brca_combined/00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz \
	    --model_db ../input/gtex_v8_eqtl_dbs_mashr/mashr_${tiss}.db \
	    --output_dir ../output/intermediate_data/gene_snp_lists/ERNEG/eqtl/${tiss}/ \
	    --output ../output/intermediate_data/gene_snp_lists/ERNEG/eqtl/${tiss}/${tiss} \
	    --parsimony 7 \
	    --gene_whitelist ../input/ERNEG_ensg_v26_ids.txt "gene_id.v26"
done
