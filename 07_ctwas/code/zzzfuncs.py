#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import parsl
import subprocess
from subprocess import call
from parsl.app.app import bash_app, python_app
from parsl.dataflow.futures import AppFuture

from zzzconfig import id_for_memo_File

@bash_app(cache=True)
def ctwas(study,
          tissue,
          eqtl_or_sqtl,
          chr_num,
          ncore,
          max_mag_snp_z,
          inputs=[],
          outputs=[],
          stdout=parsl.AUTO_LOGNAME, 
          stderr=parsl.AUTO_LOGNAME): 
    import os
    opath = outputs[0].filepath
    opath_no_gz = os.path.splitext(opath)[0]
    if os.path.exists(opath):
        return("echo 'Output exists. Remove it or delete it.'")
    filt_preharmo_z_snp = inputs[0]
    preharmo_z_snp = inputs[1]

    bash_command = \
    f"""
    Rscript 02_process_ctwas.R \
      --study {study} \
      --tissue {tissue} \
      --eqtl_or_sqtl {eqtl_or_sqtl} \
      --chr_num {chr_num} \
      --ncore {ncore} \
      --max_mag_snp_z {max_mag_snp_z} \
      --preharmo_z_snp {preharmo_z_snp} \
      --filt_preharmo_z_snp {filt_preharmo_z_snp}
    """
    return(bash_command)

@bash_app(cache=True)
def preharmonize_study(study, max_mag_snp_z,
          outputs=[],
          stdout=parsl.AUTO_LOGNAME, 
          stderr=parsl.AUTO_LOGNAME): 
    import os
    opath = outputs[0].filepath
    opath_no_gz = os.path.splitext(opath)[0]
    if os.path.exists(opath):
        return("echo 'Output exists. Remove it or delete it.'")

    bash_command = \
    f"""
    Rscript 01_preharmonize_study.R \
      --study {study} \
      --max_mag_snp_z {max_mag_snp_z} 
    """
    return(bash_command)

@bash_app(cache=True)
def combine_res(study, eqtl_or_sqtl, outputs=[], stdout=parsl.AUTO_LOGNAME,
        stderr=parsl.AUTO_LOGNAME): 
    import os
    opath = outputs[0].filepath
    opath_no_gz = os.path.splitext(opath)[0]
    if os.path.exists(opath):
        return("echo 'Output exists. Remove it or delete it.'")

    bash_command = \
    f"""
    echo "Make this file, nerd."
    """
    return(bash_command)
