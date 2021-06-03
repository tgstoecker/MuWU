rule cutadapt:
    input:
        ["RawReads/{sample}_1.fq.gz", "RawReads/{sample}_2.fq.gz"]
    output:
        fastq1="cut_reads/{sample}.fq1.gz",
        fastq2="cut_reads/{sample}.fq2.gz",
        qc="cut_reads/{sample}.qc.txt"
    params:
        adapters = "-g ^AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
                " -g ^GCCTTGGCAGTCTCAG"
                " -a GAGATAATTGCCATTATRGAMGAAGAGVG"
                " -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
                " -a ATCTCGTATGCCGTCTTCTGCTTG"
                " -G ^CAAGCAGAAGACGGCATACGAGAT"
                " -G ^GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTC"
                " -G ^CBCTCTTCKTCYATAATGGCAATTATCTC"
                " -A CTGAGACTGCCAAGGC"
                " -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        extra = "-n 8 --minimum-length 12 -e 0.2"
    log:
        "logs/cutadapt/{sample}.log"
    threads: config["threads_cutadapt"]
#    conda: "identification.yaml"
    wrapper:
        "0.74.0/bio/cutadapt/pe"


rule trimmomatic:
    input:
        r1="cut_reads/{sample}.fq1.gz",
        r2="cut_reads/{sample}.fq2.gz"
    output:
        r1="trimmed_reads/{sample}.forward_paired.fq.gz",
        r2="trimmed_reads/{sample}.reverse_paired.fq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed_reads/{sample}.forward_unpaired.fq.gz",
        r2_unpaired="trimmed_reads/{sample}.reverse_unpaired.fq.gz"
    log:
        "logs/trimmomatic/{sample}.overall.log"
    params:
        trimmer=["SLIDINGWINDOW:4:15 MINLEN:12"],
        compression_level="-9"
    threads: config["threads_trimmomatic"]
#    conda: "identification.yaml"
    wrapper:
        "0.42.0/bio/trimmomatic/pe"
