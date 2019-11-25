import re
import pandas as pd
from snakemake.utils import validate

GATK = config["gatk_cmd"]
PICARD_MAX_RECORDS = f'--MAX_RECORDS_IN_RAM {config["max_records"]}'
PICARD_TMP_DIR = f'--TMP_DIR {config["tmp_dir"]}'

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, "samples.schema.yaml")


def get_rg(wildcards):
    prefix = samples.loc[wildcards.sample]['prefix']
    match = re.match(r'[^A-Z]*([A-Z]+.[0-9]+)', prefix)
    return match.group(1)


def get_fastqs_for_sample_id(wildcards):
    prefix = samples.loc[wildcards.sample]['prefix']
    fq_dir = config['fastq_dir']
    fastqs = {'fq1': os.path.join(fq_dir, f'{prefix}.unmapped.1.fastq.gz'),
              'fq2': os.path.join(fq_dir, f'{prefix}.unmapped.2.fastq.gz')}
    return fastqs


rule all:
    input:
        expand("processed_bams/{sample}.bam", sample=samples.index)


rule bwa_map:
    input:
        config['reference'],
        unpack(get_fastqs_for_sample_id)
    output:
        temp("mapped_reads/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} 2> {log} | samtools view -b - > {output}"


rule fastq_to_ubam:
    input:
        unpack(get_fastqs_for_sample_id)
    output:
        temp("ubams/{sample}.bam")
    params:
        rg=get_rg,
        platform="illumina"
    log:
        "logs/gatk/FastqToSam/{sample}.log"
    group: "postprocessing"
    threads: 1
    shell:
        "{GATK} FastqToSam {PICARD_MAX_RECORDS} {PICARD_TMP_DIR} "
        "-F1 {input.fq1} -F2 {input.fq2} -O {output} "
        "-SM {wildcards.sample} -LB {wildcards.sample} "
        "-RG {params.rg} -PU {params.rg}.{wildcards.sample} "
        "-PL {params.platform} 2>{log}"

rule merge_ubam:
    input:
        ref=config['reference'],
        ref_dict=config['reference'].rsplit(".", 1)[0] + ".dict",
        ubam="ubams/{sample}.bam",
        bam="mapped_reads/{sample}.bam"
    output:
        temp("merged_bams/{sample}.bam")
    log:
        "logs/gatk/MergeBamAlignment/{sample}.log"
    threads: 1
    shell:
        "{GATK} MergeBamAlignment {PICARD_MAX_RECORDS} {PICARD_TMP_DIR} "
        "-R {input.ref} -O {output} "
        "-UNMAPPED {input.ubam} -ALIGNED {input.bam} 2>{log}"

rule mark_duplicates:
    input:
        "merged_bams/{sample}.bam"
    output:
        bam=temp("deduped_bams/{sample}.bam"),
        txt="metrics/{sample}.dup_metrics.txt"
    params:
        so="queryname"
    log:
        "logs/gatk/MarkDuplicates/{sample}.log"
    threads: 1
    group: "postprocessing"
    shell:
        "{GATK} MarkDuplicates {PICARD_TMP_DIR} "
        "--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 "
        "-I {input} -O {output.bam} "
        "-M {output.txt} -ASO {params.so} 2>{log}"

rule sort_bam:
    input:
        "deduped_bams/{sample}.bam"
    output:
        bam=protected("processed_bams/{sample}.bam"),
        bai=protected("processed_bams/{sample}.bai")
    params:
        so="coordinate"
    log:
        "logs/gatk/SortSam/{sample}.log"
    threads: 1
    group: "postprocessing"
    shell:
        "{GATK} SortSam {PICARD_MAX_RECORDS} {PICARD_TMP_DIR}"
        " -I {input} -O {output.bam} -SO {params.so} "
        " --CREATE_INDEX 2>{log}"
