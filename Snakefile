import pandas as pd
# from snakemake.utils import validate

GATK = config["gatk_cmd"]
PICARD_MAX_RECORDS = f'--MAX_RECORDS_IN_RAM {config["max_records"]}'
PICARD_TMP_DIR = f'--TMP_DIR {config["tmp_dir"]}'
GATK_FILTERS = ("-RF MappingQualityReadFilter --minimum-mapping-quality 30 "
                "-RF OverclippedReadFilter --filter-too-short 50")

# samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
samples = (
    pd.read_table(config["samples"])
    .set_index("sample", drop=False)
    .assign(group=lambda x: x.index.map(lambda name: '_'.join(name.split("_")[:-1])))
)
groups = list(set(samples.group))

# validate(samples, "samples.schema.yaml")


def get_rg(wildcards):
    return samples.loc[wildcards.sample]['prefix']


def get_fastqs_for_sample_id(wildcards):
    prefix = get_rg(wildcards)
    return {f'fq{i}': f'data/{prefix}.unmapped.{i}.fastq.gz' for i in '12'}


def get_samples_for_group(wildcards):
    names = samples.query(f'group=="{wildcards.group}"').index
    return [f'processed_bams/{sample}.bam' for sample in names]


rule all:
    input:
        expand("svaba/{group}.svaba.refiltered.somatic.sv.vcf", group=groups)


rule counts:
    input:
        expand("read_depth/{sample}.counts.tsv", sample=samples.index),
        expand("allelic_depth/{sample}.AD.tsv", sample=samples.index)


rule metrics:
    input:
        expand("metrics/{sample}.alignment_summary_metrics",
               sample=samples.index),


rule align:
    input:
        expand("mapped_reads/{sample}.bam", sample=samples.index)


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

# rule collect_allelic_counts:
#     input:
#         ref = config['reference'],
#         intervals = config['snp_sites'],
#         bam = "processed_bams/{sample}.bam"
#     output:
#         "allelic_depth/{sample}.AD.tsv"
#     log:
#         "logs/gatk/CollectAllelicCounts/{sample}.log"
#     shell:
#         "{GATK} CollectAllelicCounts -I {input.bam} -L {input.intervals} "
#         "-R {input.ref} -O {output} 2>{log}"

rule count_reads_allelic:
    input:
        ref = config['reference'],
        intervals = config['snp_sites'],
        bam = "processed_bams/{sample}.bam"
    output:
        "allelic_depth/{sample}.AD.tsv"
    log:
        "logs/gatk/ASEReadCounter/{sample}.log"
    shell:
        "{GATK} ASEReadCounter -I {input.bam} -V {input.intervals} "
        "-O {output} 2>{log}"

# rule collect_allelic_counts:
#     input:
#         intervals = config['snp_sites']
#         bam = "processed_bams/{sample}.bam"
#     output:
#         "allelic_depth/{sample}.AD.tsv"
#     log:
#         "logs/gatk/GetPileupSummaries/{sample}.log"
#     shell:
#         "{GATK} GetPileupSummaries -I {input.bam} -L {input.intervals} "
#         "-V {input.intervals} -O {output} 2>{log}"

rule call_structural_variants:
    input:
        ref = config['reference'],
        bam = get_samples_for_group
        # simple = config['simple_repeats'],
        # germline = config['germline_svs']
    threads: 8
    params:
        bams = lambda wildcards, input: ' '.join([f"-t {b}" for b in input.bam]),
        normal = "-n " + config['normal'],
        flags = "--min-overlap 25"
    log:
        "svaba/{group}.log"
    output:
        "svaba/{group}.svaba.unfiltered.somatic.sv.vcf"
    shell:
        "svaba run -a svaba/{wildcards.group} -p {threads} "
        "-G {input.ref} {params} "
        # "-V {input.germline} -R {input.simple}"

rule filter_structural_variants:
    input:
        "svaba/{group}.svaba.unfiltered.somatic.sv.vcf"
    threads: 1
    log:
        "logs/bcftools/{group}.log"
    output:
        "svaba/{group}.svaba.refiltered.somatic.sv.vcf"
    params:
        "-i '(SPAN<0 || SPAN>150000) && N_PASS(AD>0)==1 && AD>=3'"
    shell:
        "bcftools view {input} {params} -o {output} 2>{log}"
