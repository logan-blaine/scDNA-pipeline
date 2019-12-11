import pandas as pd
import subprocess
from snakemake.utils import validate

GATK = config["gatk_cmd"]
PICARD_MAX_RECORDS = f'--MAX_RECORDS_IN_RAM {config["max_records"]}'
PICARD_TMP_DIR = f'--TMP_DIR {config["tmp_dir"]}'

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, "samples.schema.yaml")


def get_rg(wildcards):
    return samples.loc[wildcards.sample]['prefix']


def get_fastqs_for_sample_id(wildcards):
    prefix = get_rg(wildcards)
    return {f'fq{i}': f'data/{prefix}.unmapped.{i}.fastq.gz' for i in '12'}


def gatk_inputs(wildcards, inputs):
    return ' '.join('-I '+i for i in inputs)


rule all:
    input:
        expand("processed_bams/{sample}.bam", sample=samples.index),
        expand("processed_bams/{sample}.bai", sample=samples.index)


rule counts:
    input:
        expand("read_depth/{sample}.counts.tsv", sample=samples.index)


rule metrics:
    input:
        expand("metrics/{sample}.alignment_summary_metrics",
               sample=samples.index),


rule align:
    input:
        expand("mapped_reads/{sample}.bam", sample=samples.index)


rule count:
    input:
        "cnv.pon.hdf5"


rule panel_of_normals:
    input:
        expand("read_depth/{sample}.counts.tsv", sample=samples.index)
    output:
        "cnv.pon.hdf5"
    params:
        gatk_inputs,
        "--minimum-interval-median-percentile 10",
        "--maximum-zeros-in-sample-percentage 1",
        "--maximum-zeros-in-interval-percentage 20",
        "--extreme-sample-median-percentile 10"
    log:
        "logs/gatk/CreateReadCountPanelOfNormals/PoN.log"
    shell:
        "{GATK} CreateReadCountPanelOfNormals {params} -O {output} 2> {log}"

rule bwa_map:
    input:
        config['reference'],
        unpack(get_fastqs_for_sample_id)
    output:
        temp("mapped_reads/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 16
    shell:
        "bwa mem -Y -M -t {threads} {input} 2> {log} "
        " | samtools view -b - > {output}"


rule fastq_to_ubam:
    input:
        unpack(get_fastqs_for_sample_id)
    output:
        temp("ubams/{sample}.bam")
    params:
        rg = get_rg,
        platform = "illumina"
    log:
        "logs/gatk/FastqToSam/{sample}.log"
    shell:
        "{GATK} FastqToSam {PICARD_MAX_RECORDS} {PICARD_TMP_DIR} "
        "-F1 {input.fq1} -F2 {input.fq2} -O {output} "
        "-SM {wildcards.sample} -LB {wildcards.sample} "
        "-RG {params.rg} -PU {params.rg} "
        "-PL {params.platform} 2>{log}"


rule merge_ubam:
    input:
        ref = config['reference'],
        ref_dict = config['reference'].rsplit(".", 1)[0] + ".dict",
        ubam = "ubams/{sample}.bam",
        bam = "mapped_reads/{sample}.bam"
    output:
        temp("merged_bams/{sample}.bam")
    group:
        "postprocessing"
    log:
        "logs/gatk/MergeBamAlignment/{sample}.log"
    params:
        "--ALIGNER_PROPER_PAIR_FLAGS",
        "-SO unsorted",
        "-MAX_GAPS -1"
    shell:
        "{GATK} MergeBamAlignment {PICARD_MAX_RECORDS} {PICARD_TMP_DIR} "
        "-R {input.ref} -O {output} {params} "
        "-UNMAPPED {input.ubam} -ALIGNED {input.bam} 2>{log}"


rule mark_duplicates:
    input:
        "merged_bams/{sample}.bam"
    output:
        bam = temp("deduped_bams/{sample}.bam"),
        txt = "metrics/{sample}.dup_metrics.txt"
    params:
        so = "queryname",
        px_dist = 2500
    group:
        "postprocessing"
    log:
        "logs/gatk/MarkDuplicates/{sample}.log"
    shell:
        "{GATK} MarkDuplicates {PICARD_TMP_DIR} "
        "--OPTICAL_DUPLICATE_PIXEL_DISTANCE {params.px_dist} "
        "-I {input} -O {output.bam} "
        "-M {output.txt} -ASO {params.so} 2>{log}"


rule sort_bam:
    input:
        "deduped_bams/{sample}.bam"
    output:
        bam = "processed_bams/{sample}.bam",
        bai = "processed_bams/{sample}.bai"
    params:
        so = "coordinate"
    group:
        "postprocessing"
    log:
        "logs/gatk/SortSam/{sample}.log"
    shell:
        "{GATK} SortSam {PICARD_MAX_RECORDS} {PICARD_TMP_DIR} "
        " -I {input} -O {output.bam} -SO {params.so} "
        "--CREATE_INDEX 2>{log}"


rule collect_metrics:
    input:
        ref = config['reference'],
        bam = "processed_bams/{sample}.bam"
    output:
        "metrics/{sample}.alignment_summary_metrics"
    params:
        "--PROGRAM null"
        "--PROGRAM CollectAlignmentSummaryMetrics",
        "--PROGRAM CollectInsertSizeMetrics",
        "--PROGRAM CollectSequencingArtifactMetrics",
        "--PROGRAM CollectGcBiasMetrics"
    log:
        "logs/gatk/CollectMultipleMetrics/{sample}.log"
    shell:
        "{GATK} CollectMultipleMetrics {PICARD_TMP_DIR} {params} "
        "-I {input.bam} -O metrics/{wildcards.sample} -R {input.ref} 2>{log}"


rule collect_read_counts:
    input:
        intervals = config['intervals'],
        bam = "processed_bams/{sample}.bam"
    output:
        "read_depth/{sample}.counts.tsv"
    params:
        "--interval-merging-rule OVERLAPPING_ONLY",
        "--format TSV"
    log:
        "logs/gatk/CollectReadCounts/{sample}.log"
    shell:
        "{GATK} CollectReadCounts -I {input.bam} -L {input.intervals} "
        "{params} -O {output} 2>{log}"
