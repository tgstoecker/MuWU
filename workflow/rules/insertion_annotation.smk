rule Assign_Gene_and_Transcript_IDs:
    input:
        "MuSeq_table_final/SLI-MuSeq_FGS.csv"
    output:
        "MuSeq_table_final/Mu_single_GeneIds_gene_lengths_and_stock.csv",
        "MuSeq_table_final/Mu_single_TranscriptIds_transcript_lengths_and_stock.csv"
    conda: "annotation.yaml"
    shell:
        "Rscript Annotation_of_Insertions.R"
