rule remove_duplicates_picard:
    input:
         "sorted_alignments/{sample}.sorted.bam"
    output:
         bam="removed_duplicates_alignments/{sample}.dedup.bam",
         txt="removed_duplicates_alignments/{sample}.dedup.txt"
    log:
         "logs/picard/{sample}.dedup.log"
    conda: "identification.yaml"
    shell:
         "picard MarkDuplicates I={input} O={output.bam} M={output.txt} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT > {log} 2>&1"
