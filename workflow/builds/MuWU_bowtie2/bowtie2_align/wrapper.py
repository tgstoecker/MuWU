from snakemake.shell import shell
import os.path
from os import path

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


# find check for if reads are SE
if path.exists("checks/read_type_is.SE"):

    shell(
        "(bowtie2 --threads {snakemake.threads} {extra} "
        "-x {snakemake.params.index} -U {snakemake.input.r1} "
        "--rg-id {snakemake.wildcards.sample} "
        "| samtools view -Sbh -o {snakemake.output[0]} -) {log}"
    )


# find check for if reads are PE
if path.exists("checks/read_type_is.PE"):

    print("Second/Paired sample is hidden from workflow at this point- this is on purpose.")
    print("For more info see: workflow/build/MuWU_cutadapt/wrapper.py")

    shell(
        "(bowtie2 --threads {snakemake.threads} {extra} "
        "-x {snakemake.params.index} -1 {snakemake.input.r1} -2 results/trimmed_reads/{snakemake.wildcards.sample}.reverse_paired.fq.gz "
        "--rg-id {snakemake.wildcards.sample} "
        "| samtools view -Sbh -o {snakemake.output[0]} -) {log}"
    )


