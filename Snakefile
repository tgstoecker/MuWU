SAMPLES, = glob_wildcards("RawReads/{sample}_1.fq")
INDEX_COUNT_LIST = range(1, 9)
configfile: "config.yaml"

rule all:
    input:
#        expand("cut_reads/{sample}_{paired}.cut.fq", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed_reads/{sample}_{direction}_{pair}.fq", sample=SAMPLES, direction=["forward", "reverse"], pair=["paired", "unpaired"]),
#        "FGS/Zea_mays.B73_RefGen_v4.44.gtf",
        "FGS/bwa_index.amb",
        "FGS/bwa_index.ann",
        "FGS/bwa_index.bwt",
        "FGS/bwa_index.pac",
        "FGS/bwa_index.sa",
        expand("mapped_and_sorted/{sample}.sorted.bam", sample=SAMPLES),
#        expand("alignments/{sample}.sam", sample=SAMPLES),
#        expand("sorted_alignments/{sample}_sorted.bam", sample=SAMPLES),
        expand("removed_duplicates_alignments/{sample}_dedup.bam", sample=SAMPLES),
        expand("removed_duplicates_alignments/{sample}_dedup.txt", sample=SAMPLES),
        expand("removed_duplicates_alignments/{sample}_dedup.bam.bai", sample=SAMPLES),
        expand("removed_duplicates_sam/{sample}_dedup.sam", sample=SAMPLES),
##        "removed_duplicates_sam/",
        "MuSeq_table",
        "MuSeq_table_final/SLI-MuSeq_FGS.csv",
        "MuSeq_table_final/SLI-MuSeq_FGS_annotated.csv"

rule cutadapt:
    input:
        "RawReads/{sample}_1.fq", "RawReads/{sample}_2.fq"
    output:
        fq1="cut_reads/{sample}_1.cut.fq",
        fq2="cut_reads/{sample}_2.cut.fq"
    threads: 4
    log:
        "logs/cutadapt/{sample}.log"
    shell:
        "cutadapt -g ^AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
                " -g ^GCCTTGGCAGTCTCAG"
                " -a GAGATAATTGCCATTATRGAMGAAGAGVG"
                " -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
                " -a ATCTCGTATGCCGTCTTCTGCTTG"
                " -G ^CAAGCAGAAGACGGCATACGAGAT"
                " -G ^GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTC"
                " -G ^CBCTCTTCKTCYATAATGGCAATTATCTC"
                " -A CTGAGACTGCCAAGGC"
                " -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
                " --cores 4 -n 4 -o {output.fq1} -p {output.fq2} --minimum-length 12 -e 0.2 {input}"
                "> {log} "


rule trimmomatic:
    input:
        cut1="cut_reads/{sample}_1.cut.fq",
        cut2="cut_reads/{sample}_2.cut.fq"
    output:
        forward_paired="trimmed_reads/{sample}_forward_paired.fq",
        forward_unpaired="trimmed_reads/{sample}_forward_unpaired.fq",
        reverse_paired="trimmed_reads/{sample}_reverse_paired.fq",
        reverse_unpaired="trimmed_reads/{sample}_reverse_unpaired.fq"
    threads: 4
    log: 
        trimlog="logs/trimmomatic/{sample}.trimlog",
        overall="logs/trimmomatic/{sample}.overall.log"
    shell:
       "trimmomatic PE -threads 4 -trimlog {log.trimlog} {input.cut1} {input.cut2}"
       " {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired}"
       " SLIDINGWINDOW:4:15 MINLEN:12 > {log.overall} 2>&1"


rule convert_gff3_to_gtf:
    input:
       "FGS/Zea_mays.B73_RefGen_v4.44.gff3"
    output:
        "FGS/Zea_mays.B73_RefGen_v4.44.gtf"
    threads: 1
    shell:
        "gffread {input} -T -o {output}"


rule bwa_index:
    input:
        expand("FGS/{genome}.fa", genome=config["genome"])
    output:
        "FGS/bwa_index.amb",
        "FGS/bwa_index.ann",
        "FGS/bwa_index.bwt",
        "FGS/bwa_index.pac",
        "FGS/bwa_index.sa"
    threads: 8
    log:
        "logs/bwa_index/bwa_index.log"
    params:
        prefix="FGS/bwa_index",
        algorithm="bwtsw"
    wrapper:
        "0.38.0/bio/bwa/index"


rule bwa_mem_and_sorting:
    input:
        reads=["trimmed_reads/{sample}_forward_paired.fq", "trimmed_reads/{sample}_reverse_paired.fq"],
        idx="FGS/bwa_index.sa"
    output:
        "mapped_and_sorted/{sample}.sorted.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="FGS/bwa_index",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}' -M -r 1 -k 12",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-@ 4"            # Extra args for samtools/picard.
    threads: 4
    wrapper:
        "0.38.0/bio/bwa/mem"


rule remove_duplicates_picard:
    input:
         "mapped_and_sorted/{sample}.sorted.bam"
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
        "samtools index {input}"


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
        "python ./Mu_insertions.py -c 4 --both -i MuSeq_table -o MuSeq_FGS.csv"


rule Assign_Gene_and_Transcript_IDs:
    input:
        "MuSeq_table_final/SLI-MuSeq_FGS.csv"
    output:
        "MuSeq_table_final/SLI-MuSeq_FGS_annotated.csv"
    shell:
        "Rscript AssignGeneandTranscriptIDs.R"


