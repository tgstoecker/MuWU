SAMPLES, = glob_wildcards("RawReads/{sample}.1.fq")
INDEX_COUNT_LIST = range(1, 9)
configfile: "config.yaml"

rule all:
    input:
        "multiqc/multiqc.html",
        expand("sorted_alignments/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand("removed_duplicates_alignments/{sample}.dedup.bam.bai", sample=SAMPLES),
        "MuSeq_table_final/Mu_single_GeneIds_gene_lengths_and_stock.csv",
        "MuSeq_table_final/Mu_single_TranscriptIds_transcript_lengths_and_stock.csv"


rule fasta_annotation_download:
    input: "FGS/"
    output:
        fa="FGS/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",
        gtf="FGS/Zea_mays.B73_RefGen_v4.50.gtf",
        gff3="FGS/Zea_mays.B73_RefGen_v4.50.gff3",
    params: #under these locations in the config file the links to the respective files are listed
        fasta=config["fasta"],
        gtf=config["gtf"],
        gff3=config["gff3"],
    conda: 
        "download.yaml"
    shell: 
        """
        wget -O - {params.fasta} | gunzip -c > {output.fa}
        wget -O - {params.gtf} | gunzip -c > {output.gtf}
        wget -O - {params.gff3} | gunzip -c > {output.gff3}
        """

rule fastqc:
    input:
        "RawReads/{sample}.{paired}.fq"
    output:
        html="fastqc/raw/{sample}_{paired}.html",
        zip="fastqc/raw/{sample}_{paired}_fastqc.zip" # suffix _fastqc.zip necessary for multiqc to find the file
    params: "-t 2"
    threads: 2
    log:
        "logs/fastqc/raw/{sample}.{paired}.log"
    conda: "identification.yaml"
    wrapper:
        "0.42.0/bio/fastqc"


rule cutadapt:
    input:
        "RawReads/{sample}.1.fq", "RawReads/{sample}.2.fq"
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
        others = "-n 8 --minimum-length 12 -e 0.2"
    log:
        "logs/cutadapt/{sample}.log"
    threads: config["threads_cutadapt"]
    conda: "identification.yaml"
    wrapper:
        "0.42.0/bio/cutadapt/pe"
                
                
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
    conda: "identification.yaml"
    wrapper:
        "0.42.0/bio/trimmomatic/pe"


rule trimmed_fastqc:
    input:
        "trimmed_reads/{sample}.{mate}_paired.fq.gz"
    output:
        html="fastqc/trimmed/{sample}.{mate}_paired_fastqc.html",
        zip="fastqc/trimmed/{sample}.{mate}_paired_fastqc.zip" # suffix _fastqc.zip necessary for multiqc to find the file
    threads: 2
    params: "-t 2"
    log:
        "logs/fastqc/trimmed/{sample}.{mate}.log"
    conda: "identification.yaml"
    wrapper:
        "0.42.0/bio/fastqc"


rule bowtie2_index:
    input:
        one=expand("FGS/{genome}.fa", genome=config["genome"]),
        two="FGS/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",
    output:
        "FGS/bowtie2_index.1.bt2",
        "FGS/bowtie2_index.2.bt2",
        "FGS/bowtie2_index.3.bt2",
        "FGS/bowtie2_index.4.bt2",
        "FGS/bowtie2_index.rev.1.bt2",
        "FGS/bowtie2_index.rev.2.bt2",
    threads: config["threads_bowtie_index"]
    log:
        "logs/bowtie2_index/bowtie2_index.log"
    conda: "identification.yaml"
    shell:
        "bowtie2-build --threads {threads} {input.one} FGS/bowtie2_index"


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
    conda: "identification.yaml"
    wrapper:
        "0.38.0/bio/bowtie2/align"

       
rule sam_to_sorted_bam:
    input:
         "mapped/{sample}.sam"
    output:
         "sorted_alignments/{sample}.sorted.bam"
    threads: config["threads_sam_to_sorted_bam"]
    conda: "identification.yaml"
    shell:
         "samtools sort -@ {threads} -O BAM {input} -o {output}"


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
         "picard MarkDuplicates I={input} O={output.bam} M={output.txt} REMOVE_DUPLICATES=true > {log} 2>&1"


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
        

rule multiqc:
    input:
        expand("fastqc/raw/{sample}_{paired}_fastqc.zip", sample=SAMPLES, paired=[1, 2]),
        expand("cut_reads/{sample}.qc.txt", sample=SAMPLES),
        expand("logs/trimmomatic/{sample}.overall.log" ,sample=SAMPLES),
        expand("fastqc/trimmed/{sample}.{mate}_paired_fastqc.zip", sample=SAMPLES, mate=["forward", "reverse"]),
        expand("logs/bowtie2_align/{sample}.log", sample=SAMPLES),
        expand("removed_duplicates_alignments/{sample}.dedup.bam", sample=SAMPLES)
    output:
        "multiqc/multiqc.html"
    params:
        "-ip"
    log:
        "logs/multiqc.log"
    conda: "identification.yaml"
    wrapper:
        "0.42.0/bio/multiqc"


rule prepare_MuSeq_table_folder:
    input:
        # Assuming this returns a list of your samples
         contigs=expand("removed_duplicates_sam/{sample}.dedup.sam", sample=SAMPLES)
    output:
        # Don't use the trailing "/" for directories in your rules
        assembly=directory("MuSeq_table")
    run:
        os.makedirs(output.assembly)
        for contig in input.contigs:
            # Better to symlink than to copy to save some space
            # Better to make relative links, fix if you want to
            abspath = os.path.abspath(contig)
            shell("ln -s {abspath} {output.assembly}")


rule Identify_Mu_insertions:
    input:
         "MuSeq_table"
    output:
        one="MuSeq_table_final/MuSeq_FGS.csv",
        two="MuSeq_table_final/SLI-MuSeq_FGS.csv"
    threads: config["threads_identify_Mu_insertions"]
    benchmark:
        "benchmarks/Mu_insertions-benchmark.txt"
    conda: "identification.yaml"
    shell:
        "python ./Mu_insertions.py -c {threads} --both -i MuSeq_table -o MuSeq_FGS.csv"


rule Assign_Gene_and_Transcript_IDs:
    input:
        "MuSeq_table_final/SLI-MuSeq_FGS.csv"
    output:
        "MuSeq_table_final/Mu_single_GeneIds_gene_lengths_and_stock.csv",
        "MuSeq_table_final/Mu_single_TranscriptIds_transcript_lengths_and_stock.csv"
    conda: "annotation.yaml"
    shell:
        "Rscript Annotation_of_Insertions.R"
