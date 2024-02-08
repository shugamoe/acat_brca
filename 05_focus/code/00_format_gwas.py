# Formats metaxcan imputed GWAS in preparation for FOCUS
# Using code modified from Alvaro Barbeira 
# Try not to run on login node since the GWAS uncompressed in memory is a tad large for login

import os
import pandas as pd
from pathlib import Path
from glob import glob

GWAS_SOURCES = glob("../input/metaxcan_gwas_imputed/intrinsic*.txt.gz")

def main(gwas_path):
    gwas_name = Path(Path(gwas_path).stem).stem
    gwas_fname = os.path.basename(gwas_path)
    print(gwas_path)
    print(gwas_name)

    d = pd.read_table(gwas_path)

    order = ["CHR", "SNP", "BP", "A1", "A2", "Z", "P",]
    rename={
        "panel_variant_id":"SNP",
        "chromosome":"CHR",
        "position":"BP",
        "effect_allele":"A1",
        "non_effect_allele":"A2",
        "zscore":"Z",
        "pvalue":"P"
        # "sample_size":"N" # Removing sample size for now, this is put in with focus munge
    }

    print(f"converting {gwas_name}")
    d = d.assign(chromosome = d.chromosome.apply(lambda x: int(x.split("chr")[1])))
    d = d.rename(columns=rename)[order]


    print(f"saving {gwas_name}")
    oname = os.path.join("../input/focus_ready_gwas/", gwas_fname)
    d.to_csv(oname, compression="gzip", index=False, sep="\t")
    pass

for gwas in GWAS_SOURCES:
    main(gwas)
