#!/bin/bash

snakemake --dag --configfile config.yaml --config samples=test.tsv | \
    dot -Tsvg > dag.svg
