rule index_sorted_bams_with_dups:
    input:
        "sorted_alignments/{sample}.sorted.bam"
    output:
        "sorted_alignments/{sample}.sorted.bam.bai"
    threads: config["threads_index_sorted_bams_with_dups"]
    conda: "identification.yaml"
    shell:
        "samtools index -@ {threads} {input}"


rule index_sorted_bams_without_dups:
    input:
        "removed_duplicates_alignments/{sample}.dedup.bam"
    output:
        "removed_duplicates_alignments/{sample}.dedup.bam.bai"
    threads: config["threads_index_sorted_bams_without_dups"]
    conda: "identification.yaml"
    shell:
        "samtools index -@ {threads} {input}"


rule convert_bam_to_sam:
    input:
        "removed_duplicates_alignments/{sample}.dedup.bam"
    output:
        "removed_duplicates_sam/{sample}.dedup.sam"
    threads: config["threads_index_sorted_bams_without_dups"]
    conda: "identification.yaml"
    shell:
        "samtools view -@ {threads} -o {output} {input}"
