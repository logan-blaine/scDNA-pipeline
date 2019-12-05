#!/bin/bash

snakemake --dag --configfile config.yaml --config samples=samples.tsv | \
    dot -Tsvg > dag.svg
