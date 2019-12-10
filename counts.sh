#!/bin/bash
#SBATCH --mail-user=LoganJ_Blaine@dfci.harvard.edu --mail-type=ALL
#SBATCH -c 2 --mem 1G

SBATCH_OPTS="--parsable -c {cluster.cpus} --mem {cluster.mem_gb}G --output {cluster.output} -x node07"

snakemake counts --keep-going --restart-times 2 --rerun-incomplete \
    -j 50 --local-cores ${SLURM_CPUS_PER_TASK:-1} --latency-wait 60 \
    --configfile config.yaml --cluster-config cluster.yaml \
    --cluster "sbatch $SBATCH_OPTS"
