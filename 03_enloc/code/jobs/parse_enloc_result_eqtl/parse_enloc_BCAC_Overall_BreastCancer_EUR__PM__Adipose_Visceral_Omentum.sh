#!/bin/bash
#PBS -N parse_enloc_BCAC_Overall_BreastCancer_EUR__PM__Adipose_Visceral_Omentum
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l mem=4gb
#PBS -l nodes=1:ppn=1

#PBS -o logs/parse_enloc_result_eqtl/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -j oe

module load gcc/6.2.0
module load python/3.5.3


cd $PBS_O_WORKDIR 

python3 -c '
import sys
import os
import gzip

I = sys.argv[1]
O = sys.argv[2]
K = sys.argv[3]

print("Loading key")
key={}
with gzip.open(K) as f:
  f.readline()
  for l in f:
    comps = l.decode().strip().split()
    key[comps[0]] = comps[1]

print("Processing file")
with open(I) as _i:
  with gzip.open(O, "w") as _o:
    _o.write("gwas_locus\tmolecular_qtl_trait\tlocus_gwas_pip\tlocus_rcp\tlead_coloc_SNP\tlead_snp_rcp\n".encode())
    for l in _i:
        comps = l.strip().split()
        if len(comps) == 5 and len(comps[4]) > 1 and comps[4][0:2] == "-1":
            comps.insert(4,"NA")
        comps[1] = key[comps[1]]
        _l = "{}\n".format("\t".join(comps))
        _o.write(_l.encode())
print("Done")
' ../output/enloc_eqtl/BCAC_Overall_BreastCancer_EUR__PM__Adipose_Visceral_Omentum.enloc.rst \
../output/parse_enloc_result_eqtl/BCAC_Overall_BreastCancer_EUR__PM__Adipose_Visceral_Omentum.enloc.rst.gz \
../input/gencode_key_name.txt.gz