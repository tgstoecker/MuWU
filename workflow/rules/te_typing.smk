#for now this rule expects paired data and GRID design
#I can comeback sometime and implement it dynamically for every scenario

rule read_typing:
    input:
        "rawreads/{sample}_{paired}.fq.gz",
#        check="checks/read_type.check",
#        fastq_files=get_fastqs_GRID,
    output:
        touch("results/te_typing/{sample}/{paired}/{te_type}.txt")
    params:
    #noice, access to nested yaml via lambda function
    # everything stays tidy
        type_seq=lambda w: config["TE_types"][w.te_type],
#    log:
#       "logs/cutadapt/{sample}.log"
    threads: config["threads_te_typing_seqkit"]
    conda: "../envs/te_typing.yaml"
    shell:
#        "echo {params.type_seq}"
        "zcat {input} | seqkit grep --by-seq -m 0 -j {threads} -p {params.type_seq} | seqkit seq --name > {output}"


rule all_types:
    output:
        touch("results/te_typing/all.types")
    params:
        all_types=lambda w: list(config["TE_types"].keys()), 
    conda: "../envs/annotation.yaml"
    script:
        "../scripts/te_type_annotation.R"
#    shell:
#        "{params.all_types}"
