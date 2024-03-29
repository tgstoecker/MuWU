### regardless of approach GRID/GENERIC
rule remove_duplicates_picard:
    input:
        "results/mapped/{sample}.sorted.bam"
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="results/dedup/{sample}.dedup.bam",
        metrics="results/dedup/{sample}.dedup.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.dedup.log"
    params:
        extra="REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=2048
    wrapper:
        "0.74.0/bio/picard/markduplicates"
