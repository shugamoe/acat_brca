from parsl.data_provider.files import File

class SetupCOJO:
    def __init__(self):
        # What studies can we run with COJO?
        self.STUDIES = [
                        # JCM Sourced Only
                        "BCAC_Overall_BreastCancer_EUR", "BCAC_Overall_BreastCancer_EUR_expression", 
                        "meta_analysis_BCAC_UKB_ovr_brca", "meta_analysis_BCAC_UKB_ovr_brca_expression",

                        "meta_analysis_BCAC_CIMBA_erneg_brca", "meta_analysis_BCAC_CIMBA_erneg_brca_expression",
                        "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch", "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch_expression",
                        "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap", "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap_expression",

                        "BCAC_ERPOS_BreastCancer_EUR", "BCAC_ERPOS_BreastCancer_EUR_expression",
                        "BCAC_ERPOS_BreastCancer_EUR_index_ermatch", "BCAC_ERPOS_BreastCancer_EUR_index_ermatch_expression",
                        "BCAC_ERPOS_BreastCancer_EUR_index_erswap", "BCAC_ERPOS_BreastCancer_EUR_index_erswap_expression",

                        # Old
                "BCAC_ER_negative", "UKB", "BCAC_UKB_meta",
                        "BCAC_Overall_expression", "BCAC_UKB_meta_expression", "ERNEG_CIMBA_BCAC", "ERNEG_CIMBA_BCAC_expression", "BCAC_ER_positive_expression",
                        

                        # Temp testing
                        "BCAC_UKB_meta_gcta_collin_prevent_mode",

                        # Subtypes?
                        "intrinsic_subtype_1",
                        "intrinsic_subtype_1_expression",
                        "intrinsic_subtype_2",
                        "intrinsic_subtype_2_expression",
                        "intrinsic_subtype_3",
                        "intrinsic_subtype_3_expression",
                        "intrinsic_subtype_4",
                        "intrinsic_subtype_4_expression",
                        "intrinsic_subtype_5",
                        "intrinsic_subtype_5_expression",
                        ]
    
        self.GENES_TO_RUN = {
            "BCAC_Overall_BreastCancer_EUR":
                (File("../input/BCAC_Overall_BreastCancer_EUR_sqtl_acat_sig.tsv"), "group"),
             
            "BCAC_Overall_BreastCancer_EUR_expression":
                (File("../input/BCAC_Overall_BreastCancer_EUR_expression_eqtl_acat_sig.tsv"), "group"),

            "meta_analysis_BCAC_UKB_ovr_brca":
                (File("../input/meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_sig.tsv"), "group"),

            "meta_analysis_BCAC_UKB_ovr_brca_expression":
                (File("../input/meta_analysis_BCAC_UKB_ovr_brca_expression_eqtl_acat_sig.tsv"), "group"),


            "meta_analysis_BCAC_CIMBA_erneg_brca": 
                (File("../input/meta_analysis_BCAC_CIMBA_erneg_brca_sqtl_acat_sig.tsv"), "group"),

            "meta_analysis_BCAC_CIMBA_erneg_brca_expression": 
                (File("../input/meta_analysis_BCAC_CIMBA_erneg_brca_expression_eqtl_acat_sig.tsv"), "group"),

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch": 
                (File("../input/meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch_sqtl_acat_sig.tsv"), "group"),

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch_expression": 
                (File("../input/meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch_expression_eqtl_acat_sig.tsv"), "group"),

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap": 
                (File("../input/meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap_sqtl_acat_sig.tsv"), "group"),

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap_expression": 
                (File("../input/meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap_expression_eqtl_acat_sig.tsv"), "group"),

            "BCAC_ERPOS_BreastCancer_EUR": 
                (File("../input/BCAC_ERPOS_BreastCancer_EUR_sqtl_acat_sig.tsv"), "group"),

            "BCAC_ERPOS_BreastCancer_EUR_expression": 
                (File("../input/BCAC_ERPOS_BreastCancer_EUR_expression_eqtl_acat_sig.tsv"), "group"),

            "BCAC_ERPOS_BreastCancer_EUR_index_ermatch": 
                (File("../input/BCAC_ERPOS_BreastCancer_EUR_index_ermatch_sqtl_acat_sig.tsv"), "group"),

            "BCAC_ERPOS_BreastCancer_EUR_index_ermatch_expression": 
                (File("../input/BCAC_ERPOS_BreastCancer_EUR_index_ermatch_expression_eqtl_acat_sig.tsv"), "group"),

            "BCAC_ERPOS_BreastCancer_EUR_index_erswap": 
                (File("../input/BCAC_ERPOS_BreastCancer_EUR_index_erswap_sqtl_acat_sig.tsv"), "group"),

            "BCAC_ERPOS_BreastCancer_EUR_index_erswap_expression": 
                (File("../input/BCAC_ERPOS_BreastCancer_EUR_index_erswap_expression_eqtl_acat_sig.tsv"), "group"),

            "intrinsic_subtype_1": 
                (File("../input/intrinsic_subtype_1_sqtl_acat_sig.tsv"), "group"),
            "intrinsic_subtype_1_expression": 
                (File("../input/intrinsic_subtype_1_eqtl_acat_sig.tsv"), "group"),

            "intrinsic_subtype_2": 
                (File("../input/intrinsic_subtype_2_sqtl_acat_sig.tsv"), "group"),
            "intrinsic_subtype_2_expression": 
                (File("../input/intrinsic_subtype_2_eqtl_acat_sig.tsv"), "group"),

            "intrinsic_subtype_3": 
                (File("../input/intrinsic_subtype_3_sqtl_acat_sig.tsv"), "group"),
            "intrinsic_subtype_3_expression": 
                (File("../input/intrinsic_subtype_3_eqtl_acat_sig.tsv"), "group"),

            "intrinsic_subtype_4": 
                (File("../input/intrinsic_subtype_4_sqtl_acat_sig.tsv"), "group"),
            "intrinsic_subtype_4_expression": 
                (File("../input/intrinsic_subtype_4_eqtl_acat_sig.tsv"), "group"),

            "intrinsic_subtype_5": 
                (File("../input/intrinsic_subtype_5_sqtl_acat_sig.tsv"), "group"),
            "intrinsic_subtype_5_expression": 
                (File("../input/intrinsic_subtype_5_eqtl_acat_sig.tsv"), "group"),
                }
        # Map study types to the imputed GWAS results Peter or myself created
        self.GWAS_KEY = {"BCAC_Overall_BreastCancer_EUR": 
              "../input/metaxcan_gwas_imputed/BCAC_Overall_BreastCancer_EUR.txt.gz",

            "BCAC_Overall_BreastCancer_EUR_expression": 
              "../input/metaxcan_gwas_imputed/BCAC_Overall_BreastCancer_EUR.txt.gz",

            "meta_analysis_BCAC_UKB_ovr_brca":
              "../input/metaxcan_gwas_imputed/meta_analysis_BCAC_UKB_ovr_brca.txt.gz",

            "meta_analysis_BCAC_UKB_ovr_brca_expression":
              "../input/metaxcan_gwas_imputed/meta_analysis_BCAC_UKB_ovr_brca.txt.gz",

            "meta_analysis_BCAC_CIMBA_erneg_brca": 
                "../input/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz",

            "meta_analysis_BCAC_CIMBA_erneg_brca_expression": 
                "../input/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz",

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch": 
                "../input/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz",

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch_expression": 
                "../input/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz",

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap": 
                "../input/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz",

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap_expression": 
                "../input/metaxcan_gwas_imputed/meta_analysis_BCAC_CIMBA_erneg_brca.txt.gz",

            "BCAC_ERPOS_BreastCancer_EUR": 
                "../input/metaxcan_gwas_imputed/BCAC_ERPOS_BreastCancer_EUR.txt.gz",

            "BCAC_ERPOS_BreastCancer_EUR_expression": 
                "../input/metaxcan_gwas_imputed/BCAC_ERPOS_BreastCancer_EUR.txt.gz",

            "BCAC_ERPOS_BreastCancer_EUR_index_ermatch": 
                "../input/metaxcan_gwas_imputed/BCAC_ERPOS_BreastCancer_EUR.txt.gz",

            "BCAC_ERPOS_BreastCancer_EUR_index_ermatch_expression": 
                "../input/metaxcan_gwas_imputed/BCAC_ERPOS_BreastCancer_EUR.txt.gz",

            "BCAC_ERPOS_BreastCancer_EUR_index_erswap": 
                "../input/metaxcan_gwas_imputed/BCAC_ERPOS_BreastCancer_EUR.txt.gz",

            "BCAC_ERPOS_BreastCancer_EUR_index_erswap_expression": 
                "../input/metaxcan_gwas_imputed/BCAC_ERPOS_BreastCancer_EUR.txt.gz",

            "intrinsic_subtype_1": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_1.txt.gz",
            "intrinsic_subtype_1_expression": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_1.txt.gz",
            "intrinsic_subtype_2": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_2.txt.gz",
            "intrinsic_subtype_2_expression": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_2.txt.gz",
            "intrinsic_subtype_3": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_3.txt.gz",
            "intrinsic_subtype_3_expression": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_3.txt.gz",
            "intrinsic_subtype_4": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_4.txt.gz",
            "intrinsic_subtype_4_expression": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_4.txt.gz",
            "intrinsic_subtype_5": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_5.txt.gz",
            "intrinsic_subtype_5_expression": 
                "../input/metaxcan_gwas_imputed/intrinsic_subtype_5.txt.gz",
            }
        
        # Map study types to the "overall" (combined across tissue) association files # containing introns, their effects, and related info
        self.MTISS_TWAS_KEY = {"BCAC_Overall_BreastCancer_EUR": 
            "../input/acat_prep_sqtl/BCAC_Overall_BreastCancer_EUR_sqtl_acat_input.csv",

                   "BCAC_Overall_BreastCancer_EUR_expression":
                     "../input/acat_prep_sqtl/BCAC_Overall_BreastCancer_EUR_sqtl_acat_input.csv",

                   "meta_analysis_BCAC_UKB_ovr_brca":
                     "../input/acat_prep_sqtl/meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_input.csv",

                   "meta_analysis_BCAC_UKB_ovr_brca_expression":
                     "../input/acat_prep_sqtl/meta_analysis_BCAC_UKB_ovr_brca_sqtl_acat_input.csv",

            "meta_analysis_BCAC_CIMBA_erneg_brca": 
                "../input/acat_prep_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca_sqtl_acat_input.csv",

            "meta_analysis_BCAC_CIMBA_erneg_brca_expression": 
                "../input/acat_prep_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca_sqtl_acat_input.csv",

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch": 
                "../input/acat_prep_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca_sqtl_acat_input.csv",

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch_expression": 
                "../input/acat_prep_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca_sqtl_acat_input.csv",

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap": 
                "../input/acat_prep_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca_sqtl_acat_input.csv",

            "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap_expression": 
                "../input/acat_prep_sqtl/meta_analysis_BCAC_CIMBA_erneg_brca_sqtl_acat_input.csv",

            "BCAC_ERPOS_BreastCancer_EUR": 
                "../input/acat_prep_sqtl/BCAC_ERPOS_BreastCancer_EUR_sqtl_acat_input.csv",

            "BCAC_ERPOS_BreastCancer_EUR_expression": 
                "../input/acat_prep_sqtl/BCAC_ERPOS_BreastCancer_EUR_sqtl_acat_input.csv",

            "BCAC_ERPOS_BreastCancer_EUR_index_ermatch": 
                "../input/acat_prep_sqtl/BCAC_ERPOS_BreastCancer_EUR_sqtl_acat_input.csv",

            "BCAC_ERPOS_BreastCancer_EUR_index_ermatch_expression": 
                "../input/acat_prep_sqtl/BCAC_ERPOS_BreastCancer_EUR_sqtl_acat_input.csv",

            "BCAC_ERPOS_BreastCancer_EUR_index_erswap": 
                "../input/acat_prep_sqtl/BCAC_ERPOS_BreastCancer_EUR_sqtl_acat_input.csv",

            "BCAC_ERPOS_BreastCancer_EUR_index_erswap_expression": 
                "../input/acat_prep_sqtl/BCAC_ERPOS_BreastCancer_EUR_sqtl_acat_input.csv",

                   "BCAC_UKB_meta":
            "/gpfs/data/gao-lab/Julian/gaolab_hub/projects/gtex_covar_extension/multi_tiss_intronxcan/input/ukb_bcac_meta_assocations.csv",
                   "BCAC_UKB_meta_gcta_collin_prevent_mode":
            "/gpfs/data/gao-lab/Julian/gaolab_hub/projects/gtex_covar_extension/multi_tiss_intronxcan/input/ukb_bcac_meta_assocations.csv",
                   "ERNEG_CIMBA_BCAC":
            "/gpfs/data/gao-lab/Julian/gaolab_hub/projects/twas_bcac_erneg/output/054/meta_analysis_intersect_BCAC_CIMBA_ERNEG_combined_assoc.csv",
                   "intrinsic_subtype_1":
                "../input/acat_prep_sqtl/intrinsic_subtype_1_sqtl_acat_input.csv",
                   "intrinsic_subtype_1_expression":
                "../input/acat_prep_sqtl/intrinsic_subtype_1_sqtl_acat_input.csv",
                   "intrinsic_subtype_2":
                "../input/acat_prep_sqtl/intrinsic_subtype_2_sqtl_acat_input.csv",
                   "intrinsic_subtype_2_expression":
                "../input/acat_prep_sqtl/intrinsic_subtype_2_sqtl_acat_input.csv",
                   "intrinsic_subtype_3":
                "../input/acat_prep_sqtl/intrinsic_subtype_3_sqtl_acat_input.csv",
                   "intrinsic_subtype_3_expression":
                "../input/acat_prep_sqtl/intrinsic_subtype_3_sqtl_acat_input.csv",
                   "intrinsic_subtype_4":
                "../input/acat_prep_sqtl/intrinsic_subtype_4_sqtl_acat_input.csv",
                   "intrinsic_subtype_4_expression":
                "../input/acat_prep_sqtl/intrinsic_subtype_4_sqtl_acat_input.csv",
                   "intrinsic_subtype_5":
                "../input/acat_prep_sqtl/intrinsic_subtype_5_sqtl_acat_input.csv",
                   "intrinsic_subtype_5_expression":
                "../input/acat_prep_sqtl/intrinsic_subtype_5_sqtl_acat_input.csv"
                }
        
        self.STUDY_SNP_LIST_PARENT = "../output/intermediate_data/study_snp_lists/"
        # Map study types to where we keep their list of gene SNPs
        self.STUDY_SNP_LIST_DIR_KEY = {study: f"{self.STUDY_SNP_LIST_PARENT}{study}/" for study in self.STUDIES}
    
        # Map study types to the imputed GWAS results
        self.OLD_SAMPLE_N_KEY = {"BCAC_Overall_BreastCancer_EUR": 228951, # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8263739/
                        "BCAC_Overall_BreastCancer_EUR_expression": 228951,
                        "meta_analysis_BCAC_UKB_ovr_brca": 195650 + 228951, # /gpfs/data/gao-lab/Julian/gaolab_hub/projects/twas_bcac_ukb/input/data/Anno_res_anyBC_full_chrom.breast.logistic.gz Has OBS_CT=195650
                        "BCAC_UKB_meta_expression_expression": 195650 + 228951,

                        "BCAC_ER_negative": 0, # TBD
                        "BCAC_ER_positive": 175475, # TBD
                        "BCAC_ER_positive_expression": 175475, # TBD
                        "UKB": 0, # TBD
                        "BCAC_UKB_meta_gcta_collin_prevent_mode": 195650 + 228951, # /gpfs/data/gao-lab/Julian/gaolab_hub/projects/twas_bcac_ukb/input/data/Anno_res_anyBC_full_chrom.breast.logistic.gz Has OBS_CT=195650
                        "BCAC_UKB_meta_expression": 195650 + 228951,
                        "ERNEG_CIMBA_BCAC": 228951 + 43000, # https://cimba.ccge.medschl.cam.ac.uk/
                        "ERNEG_CIMBA_BCAC_expression": 228951 + 43000, # https://cimba.ccge.medschl.cam.ac.uk/
                        } 

        # Map study types to the imputed GWAS results Peter created
        # effective_n = 4 / (1/n_case + 1/n_control)
        self.SAMPLE_N_KEY = {"BCAC_Overall_BreastCancer_EUR": 4 / (1/122977 + 1/105974), # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8263739/
                        "BCAC_Overall_BreastCancer_EUR_expression": 4 / (1/122977 + 1/105974),
                        "meta_analysis_BCAC_UKB_ovr_brca": 267555.6, # Guimin
                        "meta_analysis_BCAC_UKB_ovr_brca_expression": 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494), # Guimin
                        "meta_analysis_BCAC_CIMBA_erneg_brca": 4/(1/7784 + 1/7782) + 4/(1/1630 + 1/1712) + 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494),
                        "meta_analysis_BCAC_CIMBA_erneg_brca_expression": 4/(1/7784 + 1/7782) + 4/(1/1630 + 1/1712) + 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494),
                        "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch": 4/(1/7784 + 1/7782) + 4/(1/1630 + 1/1712) + 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494),
                        "meta_analysis_BCAC_CIMBA_erneg_brca_index_ermatch_expression": 4/(1/7784 + 1/7782) + 4/(1/1630 + 1/1712) + 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494),
                        "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap": 4/(1/7784 + 1/7782) + 4/(1/1630 + 1/1712) + 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494),
                        "meta_analysis_BCAC_CIMBA_erneg_brca_index_erswap_expression": 4/(1/7784 + 1/7782) + 4/(1/1630 + 1/1712) + 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494),
                        "BCAC_ERPOS_BreastCancer_EUR": 4 / (1/69501 + 1 /105974),
                        "BCAC_ERPOS_BreastCancer_EUR_expression": 4 / (1/69501 + 1 /105974),

                        "BCAC_ERPOS_BreastCancer_EUR_index_ermatch": 4 / (1/69501 + 1 /105974),
                        "BCAC_ERPOS_BreastCancer_EUR_index_ermatch_expression": 4 / (1/69501 + 1 /105974),
                        "BCAC_ERPOS_BreastCancer_EUR_index_erswap": 4 / (1/69501 + 1 /105974),
                        "BCAC_ERPOS_BreastCancer_EUR_index_erswap_expression": 4 / (1/69501 + 1 /105974),

                        "BCAC_ER_negative": 0, # TBD
                        "BCAC_ER_positive": 4 / (1/69501 + 1 /105974), # TBD
                        "BCAC_ER_positive_expression": 4 / (1/69501 + 1 /105974), # TBD
                        "UKB": 4 / (1/10534 + 1/185116), # TBD
                        "BCAC_UKB_meta_gcta_collin_prevent_mode": 267555.6, # Guimin
                        "ERNEG_CIMBA_BCAC":4/(1/7784 + 1/7782) + 4/(1/1630 + 1/1712) + 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494), # https://cimba.ccge.medschl.cam.ac.uk/
                        "ERNEG_CIMBA_BCAC_expression":4/(1/7784 + 1/7782) + 4/(1/1630 + 1/1712) + 4/(1/4480 + 1/17588) + 4/(1/7333 +1/42892) + 4/(1/9655 + 1/45494), # https://cimba.ccge.medschl.cam.ac.uk/

                        "intrinsic_subtype_1": 4 / (1 / 45253 + 1 / 91477),
                        "intrinsic_subtype_2": 4 / (1 / 6427 + 1 / 91477),
                        "intrinsic_subtype_3": 4 / (1 / 6350 + 1 / 91477),
                        "intrinsic_subtype_4": 4 / (1 / 2884 + 1 / 91477),
                        "intrinsic_subtype_5": 4 / (1 / 8602 + 1 / 91477),
                        "intrinsic_subtype_1_expression": 4 / (1 / 45253 + 1 / 91477),
                        "intrinsic_subtype_2_expression": 4 / (1 / 6427 + 1 / 91477),
                        "intrinsic_subtype_3_expression": 4 / (1 / 6350 + 1 / 91477),
                        "intrinsic_subtype_4_expression": 4 / (1 / 2884 + 1 / 91477),
                        "intrinsic_subtype_5_expression": 4 / (1 / 8602 + 1 / 91477)
                        } 
        
        # Created via: /gpfs/data/gao-lab/Julian/gaolab_hub/projects/gtex_covar_extension/multi_tiss_intronxcan/input/0_merge_phenotype_groups.py
        self.GROUPING = "../../02_acat_eqtl_sqtl/input/combine_phenotype_groups.txt.gz"
        
        # Created via /gpfs/data/gao-lab/Julian/gaolab_hub/projects/gtex_covar_extension/multi_tiss_intronxcan/input/0_merge_dbs.py
        self.MODEL_SQTL_DB = "../../02_acat_eqtl_sqtl/input/combine_sqtl.db"
        self.MODEL_EQTL_DB = "../input/combine_eqtl_for_covar_calc.db"

        self.STUDY_SNP_LIST_PARENT = "../output/intermediate_data/study_snp_lists/"
        self.COJO_MA_DIR = "../output/intermediate_data/cojo_ma_files/"
        self.COND_SNP_LIST_PARENT = "../output/intermediate_data/cond_snp_lists/"
        self.PROXY_COND_SNP_LIST_PARENT = f"{self.COND_SNP_LIST_PARENT}proxy_work/"
        self.PROXY_COND_SNP_LIST_DIR_KEY = {study: f"{self.PROXY_COND_SNP_LIST_PARENT}{study}/" for study in self.STUDIES}
        self.COND_SNP_LIST_DIR_KEY = {study: f"{self.COND_SNP_LIST_PARENT}{study}/" for study in self.STUDIES}

        self.CLUMP_COND_SNP_LIST_PARENT = "../output/intermediate_data/clump_cond_snp_lists/"
        self.CLUMP_COND_SNP_LIST_DIR_KEY = {study: f"{self.CLUMP_COND_SNP_LIST_PARENT}{study}/" for study in self.STUDIES}
    
        self.GENCODE = "../input/gencode_v26_all.txt"
        # self.REPORTED_SNP_LIST = "../input/bcac_ukb_meta_add_new_vars_GWAS_variant_annotated_James_DH_2023Jan5.tsv"
        self.REPORTED_SNP_LIST = "../input/bcac_ukb_meta_add_new_vars_GWAS_variant_annotated_James_DH_2023Jun02.tsv"
        self.COJO_OUT_PARENT = "../output/intermediate_data/cojo_output/"
        self.COJO_OUT_DIR_KEY = {study: f"{self.COJO_OUT_PARENT}{study}/" for study in self.STUDIES}

        self.COJO_IN_PARENT = "../output/intermediate_data/cojo_input/"
        self.COJO_IN_DIR_KEY = {study: f"{self.COJO_IN_PARENT}{study}/" for study in self.STUDIES}
        # 1000G EUR + GTEX V8 Variants
        # self.COJO_BFILE_PATTERN = "../input/COJO_REF_LD/jcm_hg38_panel_var_id/bfiles_merged/hg38_EUR_chr_{chr_num}_gtex_v8_1000G.phase3.genotypes.final"
        self.COJO_BFILE_PATTERN = "../input/backup_COJO_REF_LD/hg38_EUR_chr_{chr_num}_gtex_v8_1000G.phase3.genotypes.final"

# Misc. Helper functions
def filter_res(chr_num, app_fut):
    """
    Function to filter a list of DataFuture Objects. Combine with partial from functools.
    """
    import re
    has_match = re.search("chr{}_".format(chr_num), app_fut.outputs[0].filepath)
    return has_match

def filter_fp(chr_num, file_obj):
    """
    Function to filter a list of DataFuture Objects. Combine with partial from functools.
    """
    import re
    has_match = re.search("chr{}_".format(chr_num), file_obj.filepath)
    return has_match

def make_tmp_list_file(files):
    from uuid import uuid4
    import pandas as pd
    try:
        df = pd.DataFrame(data=[i.filepath for i in files])
    except AttributeError: 
        df = pd.DataFrame(data=[file_str for file_str in files])
    tmp_file = "tmp/" + str(uuid4()) + ".tmp"
    df.to_csv(tmp_file, index=False, header=False)
    return(tmp_file)

def extract_fp_strings(files):
    r_strings = []
    for fobj in files:
        r_strings.append(str(fobj.filepath))
    return(r_strings)
