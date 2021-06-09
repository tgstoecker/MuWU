#import pandas as pd

# assign directory containing sequencing reads as specified in config.yaml
#reads_dir = config["reads_dir"]




# checkpoint code to read count and specify all the outputs
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_samples(self):
        # read in generated grid sample sheet
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
#    # get list of base_name strings
        samples = samples_df["base_name"].tolist()
        return samples

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.access_samples.get(**w)

        # the magic, such as it is, happens here: we create the
        # information used to expand the pattern, using arbitrary
        # Python code.
        samples = self.get_samples()

        pattern = expand(self.pattern, sample=samples, **w)
        return pattern


#### ONLY FOR GRID DESIGN (STOCK MATRIX) BASED ANALYSIS ####
if config["approach"] == "GRID":

    rule create_grid_sample_sheet:
        output:
            "config/grid_sample_sheet.tsv",
            "config/read_type.yaml",
        script:
            "../scripts/grid_sample_preprocessing.R"


    checkpoint access_samples:
        input:
            "config/grid_sample_sheet.tsv",
        output:
            touch("checks/checkpoint.touch"),


    ##PE_samples
    def get_PE_fastqs_GRID(wildcards):
        """Get raw FASTQ files based on automatically generated grid_sample_sheet.tsv."""
        rd = reads_dir
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
        s = samples_df.loc[ (wildcards.sample), ["base_name", "fq_1_end", "fq_2_end"] ].dropna()
        return [ f"{rd}/{s.base_name}{s.fq_1_end}", f"{rd}/{s.base_name}{s.fq_2_end}" ]


    ##PE_samples
#    def get_PE_fastqs_GRID(wildcards):
#        """Get raw FASTQ files based on automatically generated grid_sample_sheet.tsv."""
#    # read in generated grid sample sheet
#        samples = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
#    # get list of base_name strings
#        SAMPLES = samples["base_name"].tolist()
#        rd = reads_dir
#        s = samples.loc[ (wildcards.sample), ["base_name", "fq_1_end", "fq_2_end"] ].dropna()
#        return [ f"{rd}/{s.base_name}{s.fq_1_end}", f"{rd}/{s.base_name}{s.fq_2_end}" ]


#### ONLY FOR GENERAL ANALYSIS - NO STOCK MATRIX ####
elif config["approach"] == "GENERAL":
    print("GENERAL")

else:
    print("ERROR MESSAGE")
