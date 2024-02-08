echo Resetting all affected ../output/intermediate_data for: ${1}

rm -rf ../output/intermediate_data/clump_cond_snp_lists/${1}
rm -rf ../output/intermediate_data/cojo_output/${1}
rm -rf ../output/intermediate_data/cojo_input/${1}_blacklist.tsv
rm -rf ../output/intermediate_data/cojo_input/${1}_whitelist.tsv
rm -rf ../output/intermediate_data/cojo_input/${1}
rm -rf ../output/intermediate_data/cojo_ma_files/${1}_gwas.ma

rm ../output/intermediate_data/cond_snp_lists/${1}_tracker.tsv
rm -rf ../output/intermediate_data/cond_snp_lists/${1}
rm -rf ../output/intermediate_data/cond_snp_lists/proxy_work/${1}
rm -rf ../output/intermediate_data/study_snp_lists/${1}
rm -rf ../output/intermediate_data/study_snp_lists/${1}_tracker.tsv
