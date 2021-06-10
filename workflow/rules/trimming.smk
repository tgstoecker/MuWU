rule test_rt:
##    input:
##        Checkpoint_GetReadType({rt}.works.txt)
##        Checkpoint_GetReadType("{rt}")
    output:
        touch("works.txt")
#        "checks/determine_read_type_check.{read_type}.txt"
    run:
        with open("config/read_type.yaml") as file:
            read_type_file = yaml.load(file, Loader=yaml.FullLoader)
            print(read_type_file["read_type"])
            rt = read_type_file["read_type"]
#            name = "test_file"
            filename = "checks/determine_read_type_check.%s.txt" % rt
            f = open(filename, "x")


#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

#    ruleorder: cutadapt_PE_GRID > cutadapt_SE_GRID


    if Checkpoint_GetReadType("{rt}") == "PE":
        rule cutadapt_PE_GRID:
            input:
                get_PE_fastqs_GRID,
            output:
                fastq1="cut_reads/{sample}.fq1.gz",
                fastq2="cut_reads/{sample}.fq2.gz",
                qc="cut_reads/{sample}.qc.txt"
            params:
                adapters = config["adapters"],
                extra = config["cutadapt_extra"],
            log:
                "logs/cutadapt/{sample}.log"
            threads: config["threads_cutadapt"]
#    conda: "identification.yaml"
            wrapper:
                "0.74.0/bio/cutadapt/pe"

#        rule final_schminal_PE:
#            input:
#                f1="cut_reads/{sample}.fq1.gz",
#                f2="cut_reads/{sample}.fq2.gz"
#            output:
#                touch("checks/{sample}.txt")



#    if Checkpoint_GetReadType("{read_type}") == "SE":
#    if path.exists("checks/determine_read_type_check.SE.txt"):
    if 1 > 0:
        print("passes check")
        rule cutadapt_SE_GRID:
            input:
                check="works.txt",
                fastq_files=get_SE_fastqs_GRID,
            output:
                fastq="cut_reads/{sample}.fq1.gz",
                qc="cut_reads/{sample}.qc.txt"
            params:
                adapters = config["adapters"],
                extra = config["cutadapt_extra"],
            log:
               "logs/cutadapt/{sample}.log"
            threads: config["threads_cutadapt"]
#    conda: "identification.yaml"
            wrapper:
                "file:workflow/builds/my_cutadapt"
#            run:
#                if path.exists("checks/determine_read_type_check.SE.txt"):
#                    print("ha")
#                    from snakemake.shell import shell
#                    #from os import path
#                    n = len(snakemake.input.fastq_files)
#                    assert n == 1, "Input must contain 1 (single-end) element."
# 
#                    extra = snakemake.params.get("extra", "")
#                    adapters = snakemake.params.get("adapters", "")
#                    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

#                    assert (
#                        extra != "" or adapters != ""
#                    ), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

#                    shell(
#                        "cutadapt"
#                        " {snakemake.params.adapters}"
#                        " {snakemake.params.extra}"
#                        " -j {snakemake.threads}"
#                        " -o {snakemake.output.fastq}"
#                        " {snakemake.input.fastq_files}"
#                        " > {snakemake.output.qc} {log}"
#                    )
#                else:
#                    print("bah")
#        if layouts[wildcards.sample] == "paired":
#            shell("paired-end command")
    else:
        print("meh")


#rule final_schminal_SE:
#    input:
#        "cut_reads/{sample}.fq1.gz"
#    output:
#         touch("checks/{sample}.txt")

#rule trimmomatic:
#    input:
#        r1="cut_reads/{sample}.fq1.gz",
#        r2="cut_reads/{sample}.fq2.gz"
#    output:
#        r1="trimmed_reads/{sample}.forward_paired.fq.gz",
#        r2="trimmed_reads/{sample}.reverse_paired.fq.gz",
#        # reads where trimming entirely removed the mate
#        r1_unpaired="trimmed_reads/{sample}.forward_unpaired.fq.gz",
#        r2_unpaired="trimmed_reads/{sample}.reverse_unpaired.fq.gz"
#    log:
#        "logs/trimmomatic/{sample}.overall.log"
#    params:
#        trimmer=["SLIDINGWINDOW:4:15 MINLEN:12"],
#        compression_level="-9"
#    threads: config["threads_trimmomatic"]
#    conda: "identification.yaml"
#    wrapper:
#        "0.42.0/bio/trimmomatic/pe"
