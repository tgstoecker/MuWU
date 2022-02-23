#for now this rule expects gzipped, paired data and GRID design
#I can comeback sometime and implement it dynamically for every scenario

rule read_typing:
    input:
        "rawreads/{sample}_{paired}.fq.gz",
#        check="checks/read_type.check",
#        fastq_files=get_fastqs_GRID,
    output:
        "results/te_typing/pre_sorting/{sample}/{paired}/{te_type}.txt",
    params:
    #noice, access to nested yaml via lambda function
    # everything stays tidy
        type_seq=lambda w: config["TE_types"][w.te_type],
        type=lambda w: config["TE_types"],
#    log:
#       "logs/cutadapt/{sample}.log"
    threads: config["threads_te_typing_seqkit"]
    conda: "../envs/te_typing.yaml"
    shell:
#        "echo {params.type_seq}"
        "zcat {input} | "
        "seqkit grep --by-seq -m 0 -j {threads} -p {params.type_seq} | " 
        "seqkit seq --name | " 
        "sed 's/\s.*$//' | "
        """awk '{{print $1,"{wildcards.paired}","{wildcards.te_type}"}}' OFS="\\t" > {output}"""


rule merging_te_typing:
    input:
        expand("results/te_typing/pre_sorting/{sample}/{paired}/{te_type}.txt", sample=SAMPLES, paired=PAIRED, te_type=TE_TYPES)
#    params:
#        type_seq=lambda w: config["TE_types"][w.te_type],
#        type=lambda w: config["TE_types"],
    output:
        "results/te_typing/pre_sorting/{sample}/{sample}_te_types_merged.tsv"
    shell:
         """
         cat results/te_typing/pre_sorting/{wildcards.sample}/1/* results/te_typing/pre_sorting/{wildcards.sample}/2/* > {output} &&
         sed  -i '1i Name\\tStrand\\tType' {output}
         """

#rule te_typing_annotation:
#    output:
#        touch("results/te_typing/all.types")
#    params:
#        samples=lambda w: list(config["SAMPLES"]),
#        all_types=lambda w: list(config["TE_types"].keys()), 
#    conda: "../envs/annotation.yaml"
#    script:
#        "../scripts/te_type_annotation.R"
#    shell:
#        "{params.all_types}"



###### once done with annotation  ###########
#extract all reads from fqs that correspond to uncategorized insertions

rule get_uncategorized_reads:
    input:
        reads="rawreads/{sample}_{paired}.fq.gz",
        unc_headers="headers_strand_{paired}_uncategorized_ins.csv",
#        check="checks/read_type.check",
#        fastq_files=get_fastqs_GRID,
    output:
        "results/te_typing/uncategorized/{sample}/{paired}/unc.fa",
#    params:
    #noice, access to nested yaml via lambda function
    # everything stays tidy
#        type_seq=lambda w: config["TE_types"][w.te_type],
#        type=lambda w: config["TE_types"],
#    log:
#       "logs/cutadapt/{sample}.log"
    threads: config["threads_te_typing_seqkit"]
    conda: "../envs/te_typing.yaml"
    shell:
        """
        zcat {input.reads} |
        seqkit grep -j {threads} --by-name -r -f {input.unc_headers} | 
        seqkit grep -j {threads} --by-seq -m 4 -p ATAATGGCAATTATCTC |
        seqkit seq --seq-type DNA --reverse --complement |
        seqkit fq2fa > {output}
        """
