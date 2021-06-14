#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

    rule create_read_type_check:
        output:
            touch("checks/read_type.check")
        run:
            with open("config/read_type.yaml") as file:
                read_type_file = yaml.load(file, Loader=yaml.FullLoader)
                rt = read_type_file["read_type"]
                filename = "checks/read_type_is.%s" % rt
                f = open(filename, "x")


    rule cutadapt_GRID:
        input:
            check="checks/read_type.check",
            fastq_files=get_fastqs_GRID,
        output:
            fastq1="results/cut_reads/{sample}.fq1.gz",
            # if samples are paired-end the second file is referred to in the wrapper
            qc="results/cut_reads/{sample}.qc.txt"
        params:
            adapters = config["adapters"],
            extra = config["cutadapt_extra"],
        log:
           "logs/cutadapt/{sample}.log"
        threads: config["threads_cutadapt"]
#    conda: "identification.yaml"
        wrapper:
            "file:workflow/builds/MuWU_cutadapt"


    rule trimmomatic_GRID:
        input:
            check="checks/read_type.check",
            r1="results/cut_reads/{sample}.fq1.gz",
#        r2="cut_reads/{sample}.fq2.gz"
        output:
            r1="results/trimmed_reads/{sample}.forward_paired.fq.gz",
#            r2="trimmed_reads/{sample}.reverse_paired.fq.gz",
            # reads where trimming entirely removed the mate
#            r1_unpaired="trimmed_reads/{sample}.forward_unpaired.fq.gz",
#            r2_unpaired="trimmed_reads/{sample}.reverse_unpaired.fq.gz"
        log:
            "logs/trimmomatic/{sample}.overall.log"
        params:
            trimmer=["SLIDINGWINDOW:4:15 MINLEN:12"],
            compression_level="-9"
        threads: config["threads_trimmomatic"]
        wrapper:
            "file:workflow/builds/MuWU_trimmomatic"
