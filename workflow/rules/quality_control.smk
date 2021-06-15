#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

    rule fastqc_GRID:
        input:
            fastq_files_1 = get_fastqs_fastqc_1_GRID,
            fastq_files_2 = get_fastqs_fastqc_2_GRID,
            check = "checks/read_type.check",
        output:
            html="results/fastqc/raw/{sample}_1_fastqc.html",
            zip="results/fastqc/raw/{sample}_1_fastqc.zip" # suffix _fastqc.zip necessary for multiqc to find the file
        params: "--quiet"
        threads: config["threads_fastqc"]
        log:
            "logs/fastqc/raw/{sample}.log"
        wrapper:
            "file:workflow/builds/MuWU_fastqc"


    rule trimmed_fastqc_GRID:
        input:
            fastq_files_1 = "results/trimmed_reads/{sample}.forward_paired.fq.gz",
        output:
            html="results/fastqc/trimmed/{sample}_forward_paired_fastqc.html",
            zip="results/fastqc/trimmed/{sample}_forward_paired_fastqc.zip" # suffix _fastqc.zip necessary for multiqc to find the file
        params: "--quiet"
        threads: config["threads_fastqc"]
        log:
            "logs/fastqc/trimmed/{sample}.log"
        wrapper:
            "file:workflow/builds/MuWU_trimmed_fastqc"


    rule multiqc_GRID:
        input:
            Checkpoint_ReadSampleSheet_GRID("results/fastqc/raw/{sample}_1_fastqc.zip"),
            Checkpoint_ReadSampleSheet_GRID("results/fastqc/trimmed/{sample}_forward_paired_fastqc.zip"),
            Checkpoint_ReadSampleSheet_GRID("results/cut_reads/{sample}.qc.txt"),
            Checkpoint_ReadSampleSheet_GRID("logs/trimmomatic/{sample}.overall.log"),
            Checkpoint_ReadSampleSheet_GRID("logs/bowtie2_align/{sample}.log"),
            Checkpoint_ReadSampleSheet_GRID("results/dedup/{sample}.dedup.bam"),
        output:
            "results/multiqc/multiqc.html"
        params:
            "-ip"
        log:
            "logs/multiqc/multiqc.log"
        wrapper:
            "0.74.0/bio/multiqc"
