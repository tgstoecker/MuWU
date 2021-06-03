rule prepare_MuSeq_table_folder:
    input:
        # Assuming this returns a list of your samples
         contigs=expand("removed_duplicates_sam/{sample}.dedup.sam", sample=SAMPLES)
    output:
        # Don't use the trailing "/" for directories in your rules
        assembly=directory("MuSeq_table")
    run:
        os.makedirs(output.assembly)
        for contig in input.contigs:
            # Better to symlink than to copy to save some space
            # Better to make relative links, fix if you want to
            abspath = os.path.abspath(contig)
            shell("ln -s {abspath} {output.assembly}")


rule Identify_Mu_insertions:
    input:
         "MuSeq_table"
    output:
        one="MuSeq_table_final/MuSeq_FGS.csv",
        two="MuSeq_table_final/SLI-MuSeq_FGS.csv"
    threads: config["threads_identify_Mu_insertions"]
    benchmark:
        "benchmarks/Mu_insertions-benchmark.txt"
    conda: "identification.yaml"
    shell:
        "python ./Mu_insertions.py -c {threads} --both -i MuSeq_table -o MuSeq_FGS.csv"
