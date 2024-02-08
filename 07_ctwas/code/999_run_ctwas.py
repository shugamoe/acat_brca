#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import glob
import re
import pandas as pd
sys.path.insert(1, 'parsl_stuff')

import parsl
from parsl.dataflow.futures import AppFuture
from zzzconfig import config, id_for_memo_File
from parsl.data_provider.files import File
from parsl.app.app import bash_app 
from parsl.app.errors import MissingOutputs, BashExitFailure
from functools import partial
from math import ceil
from datetime import datetime

from zzzfuncs import ctwas, \
        combine_res, \
        preharmonize_study


config.retries = 1
parsl.load(config)
parsl.set_stream_logger() # log lines to screen


# main <- function(study, model, model_source, method, ncore){
def run_ctwas(study, tissue, eqtl_or_sqtl, chr_num, max_mag_snp_z, ncore=14):
    # First, preharmonize the study
    harmo_name = f"{study}_max_mag_snp_z{max_mag_snp_z}"
    preharmo_output = File(f"../output/preharmo/{study}_max_mag_snp_z500_save_run.rds")
    filt_preharmo_output = File(f"../output/preharmo/{study}_max_mag_snp_z{max_mag_snp_z}_save_run.rds")
    filt_preharmo_res = preharmonize_study(study, max_mag_snp_z, outputs=[filt_preharmo_output])
            # stdout=f"{harmo_name}.out", stderr=f"{harmo_name}.err")
    preharmo_res = preharmonize_study(study, 500, outputs=[preharmo_output])
                        # stdout=f"no_filt_dummy_{harmo_name}.out", stderr=f"no_filt_dummy_{harmo_name}.err")

    # Run CTWAS (parameter estimation and finemapping steps divided inside the function)
    ctwas_name = f"{study}__max_mag_snp_z{max_mag_snp_z}__{tissue}__{eqtl_or_sqtl}_chr{chr_num}"
    output_file= File(f"../output/ctwas_rss/{ctwas_name}.susieIrss.txt")
    ctwas_res = ctwas(study=study,
          tissue=tissue,
          eqtl_or_sqtl=eqtl_or_sqtl,
          chr_num=chr_num,
          ncore=ncore,
          max_mag_snp_z=max_mag_snp_z,
          inputs=[filt_preharmo_res.outputs[0], preharmo_res.outputs[0]],
          outputs=[output_file])
          # stdout=f"{ctwas_name}.out",
          # stderr=f"{ctwas_name}.err"
          # )
    return()

db_types = ["eqtl", "sqtl"]
want_tissues = pd.read_table("../input/tissue_lists/11tiss.txt", header=None)[0].to_list()
non_intrin_studies = ["meta_bcac_ukb_ovr", "bcac_erp", "meta_bcac_cimba_erneg", "bcac_ovr"]
intrin_studies = ["intrinsic_subtype_1", "intrinsic_subtype_2",
                "intrinsic_subtype_3", "intrinsic_subtype_4", "intrinsic_subtype_5"]
all_studies = non_intrin_studies + intrin_studies
# z_mag_cutoffs = [15, 30, 500]
no_filt_cutoffs = [500]

# intrinsic subtypes 2-5 don't remove anything with abs(z) >= 15

# Testing
# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 21, 15)
# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 21, 20)
# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 21, 30)
# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 21, 50)

# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 22, 15)
# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 22, 20)
# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 22, 30)
# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 22, 50)
# run_ctwas(all_studies[0], want_tissues[0], "eqtl", 22, 500)

# Run everything, no cutoff
# for e_or_s in ['sqtl', 'eqtl']:
#     for study in all_studies:
#         for tiss in want_tissues:
#             for chr_num in range(1, 23):
#                 for z_mag in no_filt_cutoffs:
#                     run_ctwas(study, tiss, e_or_s, chr_num, z_mag)

# Run meta bcac/ukb Breast tissue only, variety of cutoffs
for e_or_s in ['eqtl']:
    for study in ['meta_bcac_ukb_ovr']:
        for tiss in ['Breast_Mammary_Tissue']:
            for chr_num in range(1, 23):
                for z_mag in [5, 10, 15, 20, 25, 30]:
                    run_ctwas(study, tiss, e_or_s, chr_num, z_mag)

# Run everything, no cutoff
# for e_or_s in ['sqtl', 'eqtl']:
#     for study in all_studies:
#         for tiss in want_tissues:
#             for chr_num in range(1, 23):
#                 for z_mag in no_filt_cutoffs:
#                     run_ctwas(study, tiss, e_or_s, chr_num, z_mag)
# 
# parsl.wait_for_current_tasks()
