from snakemake.shell import shell
import os.path
from os import path

if path.exists("checks/determine_read_type_check.SE.txt"):
    print("ha")
    from snakemake.shell import shell
    #from os import path
    n = len(snakemake.input.fastq_files)
    assert n == 1, "Input must contain 1 (single-end) element."

    extra = snakemake.params.get("extra", "")
    adapters = snakemake.params.get("adapters", "")
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    assert (
        extra != "" or adapters != ""
    ), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

    shell(
        "cutadapt"
        " {snakemake.params.adapters}"
        " {snakemake.params.extra}"
        " -j {snakemake.threads}"
        " -o {snakemake.output.fastq}"
        " {snakemake.input.fastq_files}"
        " > {snakemake.output.qc} {log}"
    )
else:
    print("bah")
