#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":
    rule bowtie2_align_GRID:
        input:
            check="checks/read_type.check",
            r1="results/trimmed_reads/{sample}.forward_paired.fq.gz",
            idx=multiext(
                "resources/genome",
                ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
                )
        output:
            "results/mapped/{sample}.sam"
        log:
            "logs/bowtie2_align/{sample}.log"
        params:
            index=lambda w, input: os.path.splitext(input.idx[0][0:16]),
            extra="-N 1"  # optional parameters
        threads: config["threads_bowtie_align"]
        wrapper:
            "file:workflow/builds/MuWU_bowtie2/bowtie2_align"


#### ONLY FOR GENERIC ANALYSIS - NO STOCK MATRIX ####
elif config["approach"] == "GENERIC":

    if not is_single_end_GENERIC_experiment(SAMPLES):

        rule bowtie2_PE_GENERIC:
            input:
                sample=["results/trimmed_reads/{sample}.1.fq.gz", "results/trimmed_reads/{sample}.2.fq.gz"],
                idx=multiext(
                    "resources/genome",
                    ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
                    )
            output:
                "results/mapped/{sample}.sam"
            log:
                "logs/bowtie2_align/{sample}.log"
            params:
                index="resources/genome",  # prefix of reference genome index (built with bowtie2-build)
                extra="-N 1"  # optional parameters
            threads: config["threads_bowtie_align"]
            wrapper:
                "file:workflow/builds/MuWU_bowtie2/bowtie2_align_GENERIC"


    if is_single_end_GENERIC_experiment(SAMPLES):

        rule bowtie2_SE_GENERIC:
            input:
                sample=["results/trimmed_reads/{sample}.1.fq.gz"],
                idx=multiext(
                    "resources/genome",
                    ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
                    )
            output:
                "results/mapped/{sample}.sam"
            log:
                "logs/bowtie2_align/{sample}.log"
            params:
                index="resources/genome",  # prefix of reference genome index (built with bowtie2-build)
                extra="-N 1"  # optional parameters
            threads: config["threads_bowtie_align"]
            wrapper:
                "file:workflow/builds/MuWU_bowtie2/bowtie2_align_GENERIC"



### regardless of GRID or GENERIC - perform samtools sort; SAM to sorted BAM ###
rule samtools_sort:
    input:
        "results/mapped/{sample}.sam"
    output:
        "results/mapped/{sample}.sorted.bam"
    log:
        "logs/samtools_sort/{sample}.samtools.sort.log"
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads:  # Samtools takes additional threads through its option -@
        config["threads_sam_to_sorted_bam"]     # This value - 1 will be sent to -@.
    wrapper:
        "0.74.0/bio/samtools/sort"
