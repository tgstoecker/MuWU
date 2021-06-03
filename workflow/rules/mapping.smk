rule bowtie2_align:
    input:
        sample=["trimmed_reads/{sample}.forward_paired.fq.gz", "trimmed_reads/{sample}.reverse_paired.fq.gz"],
        idx1="FGS/bowtie2_index.1.bt2",
        idx2="FGS/bowtie2_index.2.bt2",
        idx3="FGS/bowtie2_index.3.bt2",
        idx4="FGS/bowtie2_index.4.bt2",
        idx5="FGS/bowtie2_index.rev.1.bt2",
        idx6="FGS/bowtie2_index.rev.1.bt2"
    output:
        "mapped/{sample}.sam"
    log:
        "logs/bowtie2_align/{sample}.log"
    params:
        index="FGS/bowtie2_index",  # prefix of reference genome index (built with bowtie2-build)
        extra="-N 1"  # optional parameters
    threads: config["threads_bowtie_align"]
    conda: "bowtie2.yaml"
    shell:
        ""

rule sam_to_sorted_bam:
    input:
         "mapped/{sample}.sam"
    output:
         "sorted_alignments/{sample}.sorted.bam"
    threads: config["threads_sam_to_sorted_bam"]
    conda: "identification.yaml"
    shell:
         "samtools sort -@ {threads} -O BAM {input} -o {output}"
