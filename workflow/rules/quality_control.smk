rule fastqc:
    input:
        "rawreads/{sample}_{paired}.fq.gz"
    output:
        html="fastqc/raw/{sample}_{paired}.html",
        zip="fastqc/raw/{sample}_{paired}_fastqc.zip" # suffix _fastqc.zip necessary for multiqc to find the file
    params: "-t 2"
#    threads: 2
    log:
        "logs/fastqc/raw/{sample}.{paired}.log"
#    conda: "identification.yaml"
    wrapper:
        "v0.75.0/bio/fastqc"


#rule trimmed_fastqc:
#    input:
#        "trimmed_reads/{sample}.{mate}_paired.fq.gz"
#    output:
#        html="fastqc/trimmed/{sample}.{mate}_paired_fastqc.html",
#        zip="fastqc/trimmed/{sample}.{mate}_paired_fastqc.zip" # suffix _fastqc.zip necessary for multiqc to find the file
#    threads: 2
#    params: "-t 2"
#    log:
#        "logs/fastqc/trimmed/{sample}.{mate}.log"
#    conda: "identification.yaml"
#    wrapper:
#        "0.42.0/bio/fastqc"


#rule multiqc:
#    input:
#        expand("fastqc/raw/{sample}_{paired}_fastqc.zip", sample=SAMPLES, paired=[1, 2]),
#        expand("cut_reads/{sample}.qc.txt", sample=SAMPLES),
#        expand("logs/trimmomatic/{sample}.overall.log" ,sample=SAMPLES),
#        expand("fastqc/trimmed/{sample}.{mate}_paired_fastqc.zip", sample=SAMPLES, mate=["forward", "reverse"]),
#        expand("logs/bowtie2_align/{sample}.log", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLES)
#    output:
#        "multiqc/multiqc.html"
#    params:
#        "-ip"
#    log:
#        "logs/multiqc.log"
#    conda: "identification.yaml"
#    wrapper:
#        "0.42.0/bio/multiqc"
