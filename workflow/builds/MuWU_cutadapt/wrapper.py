from snakemake.shell import shell
import os.path
from os import path

#from snakemake.shell import shell

n = len(snakemake.input.fastq_files)

extra = snakemake.params.get("extra", "")
adapters = snakemake.params.get("adapters", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

assert (
    extra != "" or adapters != ""
), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

# find check for if reads are SE
if path.exists("checks/read_type_is.SE"):

    assert n == 1, "Input must contain 1 (single-end) element."

    shell(
        "cutadapt"
        " {snakemake.params.adapters}"
        " {snakemake.params.extra}"
        " -j {snakemake.threads}"
        " -o {snakemake.output.fastq1}"
        " {snakemake.input.fastq_files}"
        " > {snakemake.output.qc} {log}"
    )

# find check for if reads are PE
if path.exists("checks/read_type_is.PE"):

    assert n == 2, "Input must contain 2 (paired-end) elements."

    print("Second/Paired sample is hidden from workflow at this point - this is on purpose.")
    print("For more info see: workflow/build/MuWU_cutadapt/wrapper.py")

    shell(
        "cutadapt"
        " {snakemake.params.adapters}"
        " {snakemake.params.extra}"
        " -o {snakemake.output.fastq1}"
#        " -p {snakemake.output.fastq2}"
        " -p results/cut_reads/{snakemake.wildcards.sample}.fq2.gz"
        " -j {snakemake.threads}"
        " {snakemake.input.fastq_files}"
        " > {snakemake.output.qc} {log}"
    )

