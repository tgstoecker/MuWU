#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

    rule prepare_MuSeq_table_folder_GRID:
        input:
            # Assuming this returns a list of your samples
             contigs=Checkpoint_ReadSampleSheet_GRID("results/dedup_sam/{sample}.dedup.sam")
        output:
            # Don't use the trailing "/" for directories in your rules
            assembly=directory("results/MuSeq_table")
        run:
            os.makedirs(output.assembly)
            for contig in input.contigs:
                # Better to symlink than to copy to save some space
                # Better to make relative links
                abspath = os.path.abspath(contig)
                shell("ln -s {abspath} {output.assembly}")


    rule Identify_Mu_insertions:
        input:
            "results/MuSeq_table"
        output:
            one="results/MuSeq_table_final/MuSeq_FGS.csv",
            two="results/MuSeq_table_final/SLI-MuSeq_FGS.csv"
        threads: config["threads_identify_Mu_insertions"]
        benchmark:
            "benchmarks/Mu_insertions-benchmark.txt"
        conda: "../envs/identification.yaml"
        shell:
            "python workflow/scripts/Mu_insertions.py -c {threads} --both -i results/MuSeq_table -o MuSeq_FGS.csv"
