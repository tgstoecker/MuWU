# checkpoint code to assess samples in generated sheet and return them AFTER the checkpoint
class Checkpoint_ReadSampleSheet_GRID:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_samples(self):
        # read in generated grid sample sheet
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
       # get list of base_name strings
        samples = samples_df["base_name"].tolist()
        return samples
   
    def __call__(self, w):
        global checkpoints

        # wait for the results of 'config/grid_sample_sheet.tsv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.access_samples.get(**w)

        # Titus Brown rrightfully refers to this as magic;
        # information used to expand the pattern, is created using arbitrary Python code.
        samples = self.get_samples()
        pattern = expand(self.pattern, sample=samples, **w)
        return pattern



#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

    rule create_grid_sample_sheet:
        output:
            "config/grid_sample_sheet.tsv",
            "config/read_type.yaml",
        conda: 
            "../envs/gridsamplesheet.yaml"
        script:
            "../scripts/grid_sample_preprocessing.R"


    checkpoint access_samples:
        input:
            "config/grid_sample_sheet.tsv",
            "config/read_type.yaml"
        output:
            touch("checks/checkpoint.check"),



    ## call PE or SE reads depending on generated samples sheet in the GRID approach
    def get_fastqs_GRID(wildcards):
        """Get raw FASTQ files based on automatically generated grid_sample_sheet.tsv."""
        """Check whether or not single-end or paired-end reads."""
        rd = reads_dir
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
        
        # if fq_2_end is a column we continue in paired-end mode
        if 'fq_2_end' in samples_df.columns:
            s = samples_df.loc[ (wildcards.sample), ["base_name", "fq_1_end", "fq_2_end"] ].dropna()
            return [ f"{rd}/{s.base_name}{s.fq_1_end}", f"{rd}/{s.base_name}{s.fq_2_end}" ]

        # if fq_2_end is not a column we continue in single-end mode
        else:
            s = samples_df.loc[ (wildcards.sample), ["base_name", "fq_1_end"] ].dropna()
            return [ f"{rd}/{s.base_name}{s.fq_1_end}" ]    
         

        



#### ONLY FOR GENERAL ANALYSIS - NO STOCK MATRIX ####
elif config["approach"] == "GENERAL":
    print("GENERAL")

else:
    print("ERROR MESSAGE")
