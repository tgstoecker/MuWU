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
            trimmer=config["trimmer"],
            compression_level=config["compression_level"],
        threads: config["threads_trimmomatic"]
        wrapper:
            "file:workflow/builds/MuWU_trimmomatic"



#### ONLY FOR GENERIC ANALYSIS - NO STOCK MATRIX ####
elif config["approach"] == "GENERIC":

    if not is_single_end_GENERIC_experiment(SAMPLES):

        rule cutadapt_PE_GENERIC:
            input:
                get_fastqs_GENERIC,
            output:
                fastq1="results/cut_reads/{sample}.1.fq.gz",
                fastq2="results/cut_reads/{sample}.2.fq.gz",
                qc="results/cut_reads/{sample}.qc.txt"
            params:
                adapters = config["adapters"],
                extra = config["cutadapt_extra"],
            log:
                "logs/cutadapt/{sample}.log"
            threads: config["threads_cutadapt"]
            wrapper:
                "0.74.0/bio/cutadapt/pe"


        rule trimmomatic_PE_GENERIC:
            input:
                r1="results/cut_reads/{sample}.1.fq.gz",
                r2="results/cut_reads/{sample}.2.fq.gz"
            output:
                r1="results/trimmed_reads/{sample}.1.fq.gz",
                r2="results/trimmed_reads/{sample}.2.fq.gz",
            # reads where trimming entirely removed the mate
                r1_unpaired="results/trimmed_reads/{sample}.1.unpaired.fastq.gz",
                r2_unpaired="results/trimmed_reads/{sample}.2.unpaired.fastq.gz"
            log:
                "logs/trimmomatic/{sample}.overall.log"
            params:
                trimmer= config["trimmer"],
                extra="",
                compression_level=config["compression_level"],
            threads:
                config["threads_trimmomatic"]
            resources:
                mem_mb=2048
            wrapper:
                "0.74.0/bio/trimmomatic/pe"



    if is_single_end_GENERIC_experiment(SAMPLES):

        rule cutadapt_SE_GENERIC:
            input:
                get_fastqs_GENERIC,
            output:
                fastq="results/cut_reads/{sample}.1.fq.gz",
                qc="results/cut_reads/{sample}.qc.txt"
            params:
                adapters = config["adapters"],
                extra = config["cutadapt_extra"],
            log:
                "logs/cutadapt/{sample}.log"
            threads: config["threads_cutadapt"]
            wrapper:
                "0.74.0/bio/cutadapt/se"


        rule trimmomatic_SE_GENERIC:
            input:
                "results/cut_reads/{sample}.1.fq.gz",
            output:
                "results/trimmed_reads/{sample}.1.fq.gz",
            log:
                "logs/trimmomatic/{sample}.overall.log"
            params:
                trimmer= config["trimmer"],
                extra="",
                compression_level=config["compression_level"],
            threads:
                config["threads_trimmomatic"]
            resources:
                mem_mb=2048
            wrapper:
                "0.74.0/bio/trimmomatic/se"
