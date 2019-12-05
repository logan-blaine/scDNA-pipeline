#!/bin/bash
#SBATCH --mail-user=LoganJ_Blaine@dfci.harvard.edu --mail-type=ALL
#SBATCH -c 2 --mem 1G

SBATCH_OPTS="--parsable -c {cluster.cpus} --mem {cluster.mem_gb}G --output {cluster.output}"

snakemake -prk --rerun-incomplete  \
    -j 80 --local-cores $SLURM_CPUS_PER_TASK --latency-wait 60 \
    --configfile config.yaml --cluster-config cluster.yaml \
    --cluster "sbatch $SBATCH_OPTS" all
