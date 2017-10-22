#!/bin/bash

source activate dropseq

snakemake \
    -kp \
    --ri \
    -j 500 \
    --cluster-config /project2/gilad/spott/Pipelines/dropseq_pipeline/cluster.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=broadwl \
        --job-name={cluster.name} \
	--output={cluster.logfile}" \
    $*
