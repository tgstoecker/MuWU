#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

    rule create_grid_sample_sheet:
        output:
            "config/grid_sample_sheet.tsv",
            "config/read_type.yaml",
        log:
            "logs/sample_sheet/grid_sample_sheet.log"
        params:
            reads_dir=config["reads_dir"],
        conda:
            "../envs/gridsamplesheet.yaml"
        script:
            "../scripts/grid_sample_preprocessing.R"


    checkpoint access_samples:
        input:
            "config/grid_sample_sheet.tsv",
            "config/read_type.yaml"
        output:
            touch("checks/checkpoint.check")

### GENERIC workflow does not contain common rules, just functions - see common_utils.smk
