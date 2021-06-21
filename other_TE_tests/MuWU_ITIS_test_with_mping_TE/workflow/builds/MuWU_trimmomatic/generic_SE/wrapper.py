from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

# Distribute available threads between trimmomatic itself and any potential pigz instances
def distribute_threads(input_file, output_file, available_threads):
    gzipped_input_files = 1 if input_file.endswith(".gz") else 0
    gzipped_output_files = 1 if output_file.endswith(".gz") else 0
    potential_threads_per_process = available_threads // (
        1 + gzipped_input_files + gzipped_output_files
    )
    if potential_threads_per_process > 0:
        # decompressing pigz creates at most 4 threads
        pigz_input_threads = (
            min(4, potential_threads_per_process) if gzipped_input_files != 0 else 0
        )
        pigz_output_threads = (
            (available_threads - pigz_input_threads * gzipped_input_files)
            // (1 + gzipped_output_files)
            if gzipped_output_files != 0
            else 0
        )
        trimmomatic_threads = (
            available_threads
            - pigz_input_threads * gzipped_input_files
            - pigz_output_threads * gzipped_output_files
        )
    else:
        # not enough threads for pigz
        pigz_input_threads = 0
        pigz_output_threads = 0
        trimmomatic_threads = available_threads
    return trimmomatic_threads, pigz_input_threads, pigz_output_threads


def compose_input_gz(filename, threads):
    if filename.endswith(".gz") and threads > 0:
        return "<(pigz -p {threads} --decompress --stdout {filename})".format(
            threads=threads, filename=filename
        )
    return filename


def compose_output_gz(filename, threads, compression_level):
    if filename.endswith(".gz") and threads > 0:
        return ">(pigz -p {threads} {compression_level} > {filename})".format(
            threads=threads, compression_level=compression_level, filename=filename
        )
    return filename


extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)
#log = snakemake.log_fmt_shell(stdout=False, stderr=True)
compression_level = snakemake.params.get("compression_level", "-5")
trimmer = " ".join(snakemake.params.trimmer)

# Distribute threads
trimmomatic_threads, input_threads, output_threads = distribute_threads(
    snakemake.input[0], snakemake.output[0], snakemake.threads
)

# Collect files
input = compose_input_gz(snakemake.input[0], input_threads)
output = compose_output_gz(snakemake.output[0], output_threads, compression_level)

shell(
    "trimmomatic SE -threads {trimmomatic_threads} "
    "{java_opts} {extra} {input} {output} {trimmer} 2>{snakemake.log}"
)
