import os
import pandas
from glob import glob
import re

STUDY=os.path.basename(os.path.dirname(os.getcwd()))
OUT_DIR="../output/final_output/"

EG_FNAME="eqtl_focus_11tiss_meta_analysis_BCAC_CIMBA_erneg_brca_Chr4__PM__adipose_subcutaneous.focus.tsv"


TISSUES = [
'adipose_subcutaneous',
'adipose_visceral_omentum',
'breast_mammary_tissue',
'cells_cultured_fibroblasts',
'cells_ebv-transformed_lymphocytes',
'liver',
'ovary',
'spleen',
'uterus',
'vagina',
'whole_blood'
]

def main(eqtl_or_sqtl, tiss_tag="_11tiss_"):
    if eqtl_or_sqtl == "eqtl":
        pass
    elif eqtl_or_sqtl == "sqtl":
        pass
    else:
        raise Exception(f"Variable <eqtl_or_sqtl> should be one of 'eqtl' or 'sqtl', not: {eqtl_or_sqtl}")

    F=f"../output/final_output/focus_{eqtl_or_sqtl}{tiss_tag}{STUDY}_by_tissue.focus.tsv"
    d=[]
    count = 0
    for i in range(1,23): # (1-22)
        for cur_tiss in TISSUES:
            search_name = f"../output/by_tissue/{eqtl_or_sqtl}_focus{tiss_tag}{STUDY}_Chr{i}__PM__{cur_tiss}.focus.tsv" 
            #if i == 1: print(search_name)
            files_for_chrom = glob(search_name)
            if len(files_for_chrom) == 0:
                print("{}: {} missing for chr {}".format(eqtl_or_sqtl, cur_tiss, i))
            for result_file in files_for_chrom:
                fname_only = os.path.basename(result_file)
                # extract the tissue name after "__PM__" in the file
                tissue_target = re.search("__PM__(.*).focus.tsv", fname_only)[1]
                count += 1
                #print(result_file)
                d_ = pandas.read_table(result_file)
                d_[["tissue_target"]] = tissue_target
                d.append(d_)
    print("")
    print(F)
    print(count)
    print("Expect count of {} if there are 11 tissues".format(242))
    final = pandas.concat(d).drop_duplicates()
    final.to_csv(F, sep="\t", index=False)


if __name__ == "__main__":
   main("eqtl")
   main("sqtl")
   print("File one executed when ran directly")
else:
   print("File one imported")
