### regardless of approach GRID/GENERIC
rule index_sorted_bams_with_dups:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        "results/mapped/{sample}.sorted.bam.bai"
    log: "logs/index_bames/{sample}.indexing.log"
    params:
        "" # optional params string
    threads: config["threads_index_sorted_bams_without_dups"]
    wrapper:
        "0.75.0/bio/samtools/index"

rule index_sorted_bams_without_dups:
    input:
        "results/dedup/{sample}.dedup.bam"
    output:
        "results/dedup/{sample}.dedup.bam.bai"
    log: "logs/index_bames/{sample}.dedup.indexing.log"
    params:
        "" # optional params string
    threads:
        config["threads_index_sorted_bams_without_dups"]
    wrapper:
        "0.75.0/bio/samtools/index"

rule convert_bam_to_sam:
    input:
        "results/dedup/{sample}.dedup.bam"
    output:
        "results/dedup_sam/{sample}.dedup.sam"
    log:
        "logs/dedup_sam/{sample}.dedup.log"
    params:
        extra="" # optional params string
    wrapper:
        "0.75.0/bio/samtools/view"
