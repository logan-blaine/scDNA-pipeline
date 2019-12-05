#!/bin/bash

snakemake --dag --configfile test.yaml | dot -Tsvg > dag.svg
