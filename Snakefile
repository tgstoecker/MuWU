SAMPLES, = glob_wildcards("RawReads/{sample}_1.fq")
INDEX_COUNT_LIST = range(1, 9)

rule all:
    input:
#        expand("cut_reads/{sample}_{paired}.cut.fq", sample=SAMPLES, paired=[1, 2]),
#        expand("trimmed_reads/{sample}_{direction}_{pair}.fq", sample=SAMPLES, direction=["forward", "reverse"], pair=["paired", "unpaired"]),
#        "FGS/Zea_mays.B73_RefGen_v4.44.gtf",
#        "FGS/mays.ss",
#        "FGS/mays.exons",
#        expand("FGS/mays_tran.{indexCount}.ht2", indexCount=INDEX_COUNT_LIST),
#        expand("alignments/{sample}.sam", sample=SAMPLES),
#        expand("sorted_alignments/{sample}_sorted.bam", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}_dedup.bam", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}_dedup.txt", sample=SAMPLES),
#        expand("removed_duplicates_alignments/{sample}_dedup.bam.bai", sample=SAMPLES),
        expand("removed_duplicates_sam/{sample}_dedup.sam", sample=SAMPLES),
#        "removed_duplicates_sam/",
        "MuSeq_table/SLI-MuSeq_FGS.csv",
        "MuSeq_table/SLI-MuSeq_FGS_annotated.csv"

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
    log: "logs/trimmomatic/{sample}.trimlog"
    shell:
       "trimmomatic PE -threads 4 -trimlog {log} {input.cut1} {input.cut2}"
       " {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired}"
       " SLIDINGWINDOW:4:15 MINLEN:12"


rule convert_gff3_to_gtf:
    input:
       "FGS/Zea_mays.B73_RefGen_v4.44.gff3"
    output:
        "FGS/Zea_mays.B73_RefGen_v4.44.gtf"
    shell:
        "gffread {input} -T -o {output}"


rule extract_genome_splice_sites:
    input:
        "FGS/Zea_mays.B73_RefGen_v4.44.gtf"
    output:
        "FGS/mays.ss"
    shell:
        "hisat2_extract_splice_sites.py {input} > {output}"


rule extract_genome_exons:
     input:
        "FGS/Zea_mays.B73_RefGen_v4.44.gtf"
     output:
        "FGS/mays.exons"
     shell:
        "hisat2_extract_exons.py {input} > {output}"


rule hisat2_index:
     input:
         fa="FGS/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",
         splice_sites="FGS/mays.ss",
         exons="FGS/mays.exons"
     output:
          ht2_index = expand("FGS/mays_tran.{indexCount}.ht2", indexCount=INDEX_COUNT_LIST)
     threads: 16
     log:
         "logs/hisat2-index/hisat2-index.log"
     shell:
         "hisat2-build -p 16 {input.fa} --ss {input.splice_sites} --exon {input.exons} FGS/mays_tran > {log} 2>&1"


rule hisat2_mapping:
     input:
        f_paired="trimmed_reads/{sample}_forward_paired.fq",
        f_unpaired="trimmed_reads/{sample}_forward_unpaired.fq",
        r_paired="trimmed_reads/{sample}_reverse_paired.fq",
        r_unpaired="trimmed_reads/{sample}_reverse_unpaired.fq",
        index = expand("FGS/mays_tran.{indexCount}.ht2", indexCount=INDEX_COUNT_LIST)
     output:
           "alignments/{sample}.sam"
     threads: 4
     log:
        "logs/hisat2-alignments/{sample}.alignment.log"
     shell:
         "hisat2 -p 4 -x FGS/mays_tran"
         " -1 {input.f_paired} -2 {input.r_paired}"
         " -U {input.f_unpaired},{input.r_unpaired}"
         " -S {output} > {log} 2>&1"


rule sam_to_sorted_bam:
    input:
         "alignments/{sample}.sam"
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
        "samtools index {input}"


rule convert bam_to_sam:
    input:
        "removed_duplicates_alignments/{sample}_dedup.bam"
    output:
        "removed_duplicates_sam/{sample}_dedup.sam"
    threads: 16
    shell:
        "samtools view -@ 16 -o {output} {input}"


rule Identify_Mu_insertions:
    input:
#         "removed_duplicates_sam/"
#        "removed_duplicates_sam/{sample}_dedup.sam"
    output:
        "MuSeq_table/SLI-MuSeq_FGS.csv"
    threads: 16
    benchmark:
        "benchmarks/Mu_insertions-benchmark.txt"
    shell:
        "python ./Mu_insertions.py -c 16 --both -i removed_duplicates_sam/ -o MuSeq_FGS.csv"


rule Assign_Gene_and_Transcript_IDs:
    input:
        "MuSeq_table/SLI-MuSeq_FGS.csv"
    output:
        "MuSeq_table/SLI-MuSeq_FGS_annotated.csv"
    shell:
        "Rscript AssignGeneandTranscriptIDs.R"


