#!/bin/bash
#SBATCH --mail-user=LoganJ_Blaine@dfci.harvard.edu --mail-type=ALL
#SBATCH -c 2 --mem 1G

CORES=${SLURM_CPUS_PER_TASK:-1}
SAMPLES="samples.tsv"

SBATCH_OPTS="--parsable -c {cluster.cpus} --mem {cluster.mem_gb}G --output {cluster.output}"

snakemake -prk --rerun-incomplete --nolock --notemp \
    -j 50 --local-cores $CORES --latency-wait 60 \
    --configfile config.yaml --config samples=$SAMPLES \
	--cluster-config cluster.yaml \
    --cluster "sbatch $SBATCH_OPTS" align
