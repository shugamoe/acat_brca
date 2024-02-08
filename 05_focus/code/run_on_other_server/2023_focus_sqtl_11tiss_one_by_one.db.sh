module load python
# conda activate /project2/guiming/Julian/software/envs/focus
source activate /project2/guiming/Julian/software/envs/focus

data_path="/project2/guiming/Julian/FOCUS"

tissues=(Adipose_Subcutaneous Adipose_Visceral_Omentum Breast_Mammary_Tissue Cells_Cultured_fibroblasts Cells_EBV-transformed_lymphocytes Liver Ovary Spleen Uterus Vagina Whole_Blood) 

n=${#tissues[@]}

for idx in `seq 0 $((n - 1))`
do
    sleep 10
    tissue=${tissues[$idx]}

    focus import sqtl_dbs/patched/mashr_${tissue}.db predixcan \
      --tissue ${tissue} \
      --name GTEx \
      --assay rnaseq \
      --predixcan-method mashr \
      --output combined_dbs/${tissue,,}_sqtl
done
