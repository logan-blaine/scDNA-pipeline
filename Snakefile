configfile: "config.yaml"


rule all:
    input:
        expand("final_bams/{sample}.bam", sample=config['samples'])

rule add_replace_read_groups:
    input:
        lambda wildcards: expand("sorted_reads/{id}.bam", id=config['samples'][wildcards.sample])
    output:
        bam="final_bams/{sample}.bam",
        bai="final_bams/{sample}.bai"
    params:
        rg=lambda wildcards: config['samples'][wildcards.sample],
        flags="--CREATE_INDEX --VALIDATION_STRINGENCY SILENT -PL illumina"
    shell:
        "gatk AddOrReplaceReadGroups -I {input} -O {output.bam} {params.flags} "
        "-ID {params.rg} -PU {params.rg} -SM {wildcards.sample} -LB {wildcards.sample}"

rule bwa_map:
    input:
        config['reference'],
        expand("data/{{sample}}.unmapped.{lane}.fastq.gz", lane=[1,2])
    output:
        temp("mapped_reads/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 4
    shell:
        "bwa mem -t {threads} {input} | samtools view -b - > {output} 2> {log}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} {input} -o {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

