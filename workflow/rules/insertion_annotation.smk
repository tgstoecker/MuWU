#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

    rule Annotate_all_insertions_GRID:
        input:
            all="results/insertions_table_final/all_identified_insertions.csv",
            annotation="resources/final_annotation_table",
        output:
            all="results/insertions_table_final/all_identified_insertions_annotated.csv",
            germinal="results/insertions_table_final/germinal_identified_insertions_annotated.csv"
        conda: "../envs/annotation.yaml"
        script:
            "workflow/scripts/annotation_all_insertions.R"

    rule Annotate_germinal_insertions_GRID:
        input:
            germinal="results/insertions_table_final/germinal_insertions.csv",
            grid_table="config/grid_sample_sheet.tsv",
            annotation="resources/final_annotation_table",
        output:
            all="results/insertions_table_final/all_identified_insertions_annotated.csv",
            germinal="results/insertions_table_final/germinal_identified_insertions_annotated.csv"
        conda: "../envs/annotation.yaml"
        script:
            "workflow/scripts/annotation_germinal_insertions.R"


#### ONLY FOR GENERIC APPROACH BASED ANALYSIS ####
if config["approach"] == "GENERIC":

    rule Annotate_all_insertions_GENERIC:
        input:
            all="results/insertions_table_final/all_identified_insertions.csv",
            annotation="resources/final_annotation_table",
        output:
            all="results/insertions_table_final/all_identified_insertions_annotated.csv",
        conda: "../envs/annotation.yaml"
        script:
            "../scripts/annotation_all_insertions.R"
