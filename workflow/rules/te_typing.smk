#for now this rule expects paired data and GRID design
#I can comeback sometime and implement it dynamically for every scenario

rule read_typing:
#    input:
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
#    threads: config["threads_cutadapt"]
#    conda: "identification.yaml"
    shell:
        "echo {params.type_seq}"
#        "seqkit grep --by-seq -m 0 -p GTCTCTTCTTCCATAATGGCAATTATCTC | seqkit seq --name > {output}"
