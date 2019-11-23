gatk = 'gatk --java-options "-Xmx4G"'


def get_rg(wildcards):
    return config['samples'][wildcards.sample]


def get_fastqs_for_sample_id(wildcards):
    rg = get_rg(wildcards)
    return {'fq' + i: f'data/{rg}.unmapped.{i}.fastq.gz' for i in ['1', '2']}


rule all:
    input:
        expand("processed_bams/{sample}.bam", sample=config['samples'])

rule fastq_to_ubam:
    input:
        unpack(get_fastqs_for_sample_id)
    output:
        "ubams/{sample}.bam"
    params:
        rg = get_rg,
        platform = "illumina"
    shell:
        "{gatk} FastqToSam -F1 {input.fq1} -F2 {input.fq2} -O {output} "
        "-SM {wildcards.sample} -LB {wildcards.sample} "
        "-RG {params.rg} -PU {params.rg} -PL {params.platform}"

rule bwa_map:
    input:
        config['reference'],
        unpack(get_fastqs_for_sample_id)
    output:
        "mapped_reads/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 4
    shell:
        "bwa mem -t {threads} {input} 2> {log} | samtools view -b - > {output}"

rule merge_ubam:
    input:
        ref = config['reference'],
        ubam = "ubams/{sample}.bam",
        bam = "mapped_reads/{sample}.bam"
    output:
        "merged_bams/{sample}.bam"
    shell:
        "{gatk} MergeBamAlignment -R {input.ref} "
        "-UNMAPPED {input.ubam} -ALIGNED {input.bam} -O {output}"

rule mark_duplicates:
    input:
        "merged_bams/{sample}.bam"
    output:
        bam = "processed_bams/{sample}.bam",
        txt = "metrics/{sample}.dup_metrics.txt"
    params:
        so = "queryname"
    shell:
        "{gatk} MarkDuplicates -I {input} -O {output.bam} "
        "-M {output.txt} -ASO {params.so}"

# rule add_replace_read_groups:
#     input:
#         lambda wildcards: expand(
#             "sorted_reads/{id}.bam", id=config['samples'][wildcards.sample])
#     output:
#         bam = "final_bams/{sample}.bam",
#         bai = "final_bams/{sample}.bai"
#     params:
#         rg=get_rg,
#         flags = "--CREATE_INDEX --VALIDATION_STRINGENCY SILENT -PL illumina"
#     shell:
#         "gatk AddOrReplaceReadGroups -I {input} -O {output.bam} {params.flags} "
#         "-ID {params.rg} -PU {params.rg} -SM {wildcards.sample} -LB {wildcards.sample}"


# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} {input} -o {output}"
#
# rule samtools_index:
#     input:
#         "sorted_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam.bai"
#     shell:
#         "samtools index {input}"

