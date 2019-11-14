SAMPLES, = glob_wildcards("RawReads/{sample}_1.fq")
INDEX_COUNT_LIST = range(1, 9)
configfile: "config.yaml"

rule all:
    input:
#        expand("cut_reads/{sample}_{paired}.cut.fq", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed_reads/{sample}_{direction}_{pair}.fq", sample=SAMPLES, direction=["forward", "reverse"], pair=["paired", "unpaired"]),
#        "FGS/Zea_mays.B73_RefGen_v4.44.gtf",
#        "FGS/bwa_index.amb",
#        "FGS/bwa_index.ann",
#        "FGS/bwa_index.bwt",
#        "FGS/bwa_index.pac",
#        "FGS/bwa_index.sa",
#        expand("mapped_and_sorted/{sample}.sorted.bam", sample=SAMPLES),
#        expand("alignments/{sample}.sam", sample=SAMPLES),
#        expand("sorted_alignments/{sample}_sorted.bam", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}_dedup.bam", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}_dedup.txt", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}_dedup.bam.bai", sample=SAMPLES),
#        expand("removed_duplicates_sam/{sample}_dedup.sam", sample=SAMPLES),
##        "removed_duplicates_sam/",
#        "MuSeq_table",
#        "MuSeq_table_final/SLI-MuSeq_FGS.csv",
        "MuSeq_table_final/SLI-MuSeq_FGS_annotated.csv"

rule cutadapt:
    input:
        "RawReads/{sample}_1.fq", "RawReads/{sample}_2.fq"
    output:
        fastq1="cut_reads/{sample}.fq1.gz",
        fastq2="cut_reads/{sample}.fq2.gz",
        qc="cut_reads/{sample}.qc.txt"
    params:
        adapters = "g ^AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
                " -g ^GCCTTGGCAGTCTCAG"
                " -a GAGATAATTGCCATTATRGAMGAAGAGVG"
                " -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
                " -a ATCTCGTATGCCGTCTTCTGCTTG"
                " -G ^CAAGCAGAAGACGGCATACGAGAT"
                " -G ^GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTC"
                " -G ^CBCTCTTCKTCYATAATGGCAATTATCTC"
                " -A CTGAGACTGCCAAGGC"
                " -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
        others = "-n 4 --minimum-length 12 -e 0.2"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 4 
    wrapper:
        "0.42.0/bio/cutadapt/pe"
                
                
rule trimmomatic_pe:
    input:
        r1="cut_reads/{sample}.fq1.gz",
        r2="cut_reads/{sample}.fq2.gz"
    output:
        r1="trimmed_reads/{sample}_forward_paired.fq.gz",
        r2="trimmed_reads/{sample}_reverse_paired.fq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed_reads/{sample}_forward_unpaired.fq.gz",
        r2_unpaired="trimmed_reads/{sample}_reverse_unpaired.fq.gz"
    log:
        trimlog="logs/trimmomatic/{sample}.trimlog",
        overall="logs/trimmomatic/{sample}.overall.log"
    params:
        # list of trimmers (see manual)
        trimmer=["SLIDINGWINDOW:4:15 MINLEN:12"],
        # optional parameters
        extra="-trimlog {log.trimlog}",
        compression_level="-9"
    threads:
        4
    wrapper:
        "0.42.0/bio/trimmomatic/pe"
                
            
rule convert_gff3_to_gtf:
    input:
       expand("FGS/{annotation}.gff3", annotation=config["annotation"])
    output:
        "FGS/Zea_mays.B73_RefGen_v4.gtf"
    threads: 1
    shell:
        "gffread {input} -T -o {output}"


rule bowtie2_index:
    input:
        expand("FGS/{genome}.fa", genome=config["genome"])
    output:
        "FGS/bowtie2_index.1.bt2",
        "FGS/bowtie2_index.2.bt2",
        "FGS/bowtie2_index.3.bt2",
        "FGS/bowtie2_index.4.bt2",
        "FGS/bowtie2_index.rev.1.bt2",
        "FGS/bowtie2_index.rev.2.bt2",
    threads: 16
    log:
        "logs/bowtie2_index/bowtie2_index.log"
    shell:
        "bowtie2-build --threads 16 {input} FGS/bowtie2_index"


rule bowtie2:
    input:
        sample=["trimmed_reads/{sample}_forward_paired.fq.gz", "trimmed_reads/{sample}_forward_paired.fq.gz"],
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
    threads: 4
    wrapper:
        "0.38.0/bio/bowtie2/align"

       
rule sam_to_sorted_bam:
    input:
         "mapped/{sample}.sam"
    output:
         "sorted_alignments/{sample}_sorted.bam"
    threads: 4
    shell:
         "samtools sort -@ 4 -O BAM {input} -o {output}"


rule remove_duplicates_picard:
    input:
         "sorted_alignments/{sample}_sorted.bam"
    output:
         bam="removed_duplicates_alignments/{sample}_dedup.bam",
         txt="removed_duplicates_alignments/{sample}_dedup.txt"
    threads: 4
    log:
         "logs/picard/{sample}_dedup.log"
    shell:
         "picard MarkDuplicates I={input} O={output.bam} M={output.txt} REMOVE_DUPLICATES=true > {log} 2>&1"


rule index_final_bam:
    input:
        "removed_duplicates_alignments/{sample}_dedup.bam"
    output:
        "removed_duplicates_alignments/{sample}_dedup.bam.bai"
    threads: 4
    shell:
        "samtools index -@ 4 {input}"


rule convert bam_to_sam:
    input:
        "removed_duplicates_alignments/{sample}_dedup.bam"
    output:
        "removed_duplicates_sam/{sample}_dedup.sam"
    threads: 4
    shell:
        "samtools view -@ 4 -o {output} {input}"


rule prepare_MuSeq_table_folder:
    input:
        # Assuming this returns a list of your samples
        #contigs="removed_duplicates_sam/{sample}_dedup.sam"
         contigs=expand("removed_duplicates_sam/{sample}_dedup.sam", sample=SAMPLES)
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
#         "removed_duplicates_sam/"
#        "removed_duplicates_sam/{sample}_dedup.sam"
         "MuSeq_table"
    output:
        one="MuSeq_table_final/MuSeq_FGS.csv",
        two="MuSeq_table_final/SLI-MuSeq_FGS.csv"
    threads: 4
    benchmark:
        "benchmarks/Mu_insertions-benchmark.txt"
    shell:
        "python ./Mu_insertions.py -c 16 --both -i MuSeq_table -o MuSeq_FGS.csv"


rule Assign_Gene_and_Transcript_IDs:
    input:
        "MuSeq_table_final/SLI-MuSeq_FGS.csv",
        "FGS/Zea_mays.B73_RefGen_v4.gtf"
    output:
        "MuSeq_table_final/SLI-MuSeq_FGS_annotated.csv"
    shell:
        "Rscript AssignGeneandTranscriptIDs.R"
