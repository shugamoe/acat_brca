# for study_idx in {1..9}; do
#   echo ${study_idx}
  # 1=eqtl, 2=sqtl
#   for qtl_idx in {1..2}; do
#     echo ${qtl_idx}
#     sbatch --job-name=brca_coloc_${study_idx}_${qtl_idx} 999_run_coloc.sbatch ${study_idx} ${qtl_idx} 
#   done
#   sleep 1.5h;
# done

sbatch --job-name=brca_coloc_4_1 999_run_coloc.sbatch 4 1
sleep 5s;
sbatch --job-name=brca_coloc_5_2 999_run_coloc.sbatch 5 2
sleep 5s;
sbatch --job-name=brca_coloc_6_2 999_run_coloc.sbatch 6 2
sleep 5s;
sbatch --job-name=brca_coloc_6_2 999_run_coloc.sbatch 7 2
