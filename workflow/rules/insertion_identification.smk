#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

    rule prepare_MuSeq_table_folder_GRID:
        input:
            # Assuming this returns a list of your samples
             contigs=Checkpoint_ReadSampleSheet_GRID("results/dedup_sam/{sample}.dedup.sam")
        output:
            # Don't use the trailing "/" for directories in your rules
            assembly=directory("results/insertions_table")
        run:
            os.makedirs(output.assembly)
            for contig in input.contigs:
                # Better to symlink than to copy to save some space
                # Better to make relative links
                abspath = os.path.abspath(contig)
                shell("ln -s {abspath} {output.assembly}")


    rule Identify_Mu_insertions:
        input:
            file_loc="results/insertions_table",
            grid_sample_sheet="config/grid_sample_sheet.tsv",
        output:
            one="results/insertions_table_final/all_identified_insertions.csv",
            two="results/insertions_table_final/germinal_insertions.csv",
        params:
            # overlap_size refers to bp length of 
            # in the MuSeq approach as part of BonnMu this is 9 - as is characteristic for Mutator transposons
            overlap_size = config["overlap_size"],
            overlap_support = config["overlap_support"],
            approach = config["approach"],
        threads: config["threads_identify_Mu_insertions"]
        conda: "../envs/identification.yaml"
        script:
            "../scripts/insertions.py"
