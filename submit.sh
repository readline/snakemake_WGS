#!/bin/bash
snakedir=`pwd`
module load snakemake || exit 1

rm -rf logs/slurm
mkdir -p logs/slurm
sbcmd="sbatch --cpus-per-task={threads} --mem={resources.mem}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out="
sbcmd+=$snakedir
sbcmd+="/logs/slurm/{rule}.%j.out {resources.extra}"

snakemake -pr --keep-going --local-cores 4 \
    --jobs 999 --cluster-config cluster.yaml --cluster "$sbcmd" \
    --latency-wait 120 all --configfile config.yaml --snakefile Snakefile
