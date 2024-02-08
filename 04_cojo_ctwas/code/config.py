#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from parsl.config import Config
from parsl.data_provider.files import File
from parsl.channels import LocalChannel
from parsl.executors import HighThroughputExecutor
from parsl.providers import TorqueProvider, SlurmProvider
from parsl.addresses import address_by_route, address_by_query, address_by_hostname
from parsl.launchers import SingleNodeLauncher, SimpleLauncher
from parsl.utils import get_all_checkpoints
from parsl.dataflow.memoization import id_for_memo

config = Config(
    usage_tracking = True,
    checkpoint_mode = "task_exit",
    checkpoint_files = get_all_checkpoints(),
    executors= [
        HighThroughputExecutor(
            label="gp",
            max_workers=1, 
            address=address_by_hostname(), # Gets cri16in001 as address
            worker_debug = True,
            # provider=TorqueProvider(
            provider=SlurmProvider(
                channel=LocalChannel(), 
                launcher=SingleNodeLauncher(),
                worker_init="\n".join((
                    "cd $SLURM_SUBMIT_DIR",
                    "module load gcc/12.1.0;",
                    "module load miniconda3/23.1.0;"
                    "module load plink/1.9;",
                    "source activate /gpfs/data/gao-lab/Julian/software/envs/metaxcan;",
                    "export PATH=/gpfs/data/gao-lab/Julian/software/envs/metaxcan/bin:$PATH; ",
                    "export PATH=/gpfs/data/gao-lab/Julian/software/gcta-1.94.1-linux-kernel-3-x86_64:$PATH; ",
                    )),
                exclusive=False,
                nodes_per_block=1,
                cores_per_node=1,
                walltime="72:00:00",
                mem_per_node=16,
                init_blocks=0,
                min_blocks=0,
                max_blocks=40,
                parallelism=1
            )
        ),
    ]
)


@id_for_memo.register(File)
def id_for_memo_File(f, output_ref=False):
    import os
    if output_ref:
        # logger.debug("hashing File as output ref without content: {}".format(f))
        print("hashing File as output ref without content: {}".format(f))
        return f.url
    else:
        # logger.debug("hashing File as input with content: {}".format(f))
        print("hashing File as input with content: {}".format(f))
        assert f.scheme == "file"
        filename = f.filepath
        stat_result = os.stat(filename)
        
        return [f.url, stat_result.st_size, stat_result.st_mtime]
