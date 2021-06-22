from os import path
import re
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


def basename_without_ext(file_path):
    """Returns basename of file path, without the file extension."""

    base = path.basename(file_path)
    # Remove file extension(s) (similar to the internal fastqc approach)
    base = re.sub("\\.gz$", "", base)
    base = re.sub("\\.bz2$", "", base)
    base = re.sub("\\.txt$", "", base)
    base = re.sub("\\.fastq$", "", base)
    base = re.sub("\\.fq$", "", base)
    base = re.sub("\\.sam$", "", base)
    base = re.sub("\\.bam$", "", base)

    return base


# Run fastqc, since there can be race conditions if multiple jobs
# use the same fastqc dir, we create a temp dir.
with TemporaryDirectory() as tempdir:

    # find check for if reads are SE
    if path.exists("checks/read_type_is.SE"):
        shell(
            "fastqc {snakemake.params} -t {snakemake.threads} "
            "--outdir {tempdir:q} {snakemake.input.fastq_files_1:q}"
            " {log}"
        )
        # Move outputs into proper position.
        output_base = basename_without_ext(snakemake.input[0])
        html_path = path.join(tempdir, output_base + "_fastqc.html")
        zip_path = path.join(tempdir, output_base + "_fastqc.zip")

        if snakemake.output.html != html_path:
            shell("mv {html_path:q} {snakemake.output.html:q}")

        if snakemake.output.zip != zip_path:
            shell("mv {zip_path:q} {snakemake.output.zip:q}")


    # find check for if reads are PE
    if path.exists("checks/read_type_is.PE"):
        shell(
            "fastqc {snakemake.params} -t {snakemake.threads} "
            "--outdir {tempdir:q} {snakemake.input.fastq_files_1:q}"
            " {log}"
        )
        shell(
            "fastqc {snakemake.params} -t {snakemake.threads} "
            "--outdir results/fastqc/raw/ {snakemake.input.fastq_files_2:q}"
            " {log}"
        )

        # Move outputs into proper position.
        output_base_1 = basename_without_ext(snakemake.input[0])
        html_path_1 = path.join(tempdir, output_base_1 + "_fastqc.html")
        zip_path_1 = path.join(tempdir, output_base_1 + "_fastqc.zip")

        if snakemake.output.html != html_path_1:
            shell("mv {html_path_1:q} {snakemake.output.html:q}")

        if snakemake.output.zip != zip_path_1:
            shell("mv {zip_path_1:q} {snakemake.output.zip:q}")
