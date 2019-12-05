#!/bin/bash
#SBATCH --mail-user=LoganJ_Blaine@dfci.harvard.edu --mail-type=ALL

snakemake -prk --local-cores 64 -j 100 --latency-wait 60 --configfile test.yaml \
    --cluster "sbatch --parsable -c {threads} --mem 16G --output /dev/null --error /dev/null"

