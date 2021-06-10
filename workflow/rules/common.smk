# checkpoint code to read count and specify all the outputs
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

#    def __init2__(self2, pattern2):
#        self.pattern2 = pattern2

    def get_samples(self):
        # read in generated grid sample sheet
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
       # get list of base_name strings
        samples = samples_df["base_name"].tolist()
        return samples

    
#    def return_read_type(self2):
#        with open("config/read_type.yaml") as file:
#            read_type_file = yaml.load(file, Loader=yaml.FullLoader)
#            print(read_type_file["read_type"])
#            return read_type_file["read_type"]


    def __call__(self, w):
        global checkpoints

        # wait for the results of 'config/grid_sample_sheet.tsv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.access_samples.get(**w)

        # Titus Brown rrightfully refers to this as magic;
        # information used to expand the pattern, is created using arbitrary Python code.
        samples = self.get_samples()

        pattern = expand(self.pattern, sample=samples, **w)

#        rt = self2.return_read_type()
#        pattern2 = rt
      
        return pattern


#    def __call2__(self2, w):
#        global checkpoints

        # wait for the results of 'config/grid_sample_sheet.tsv'; this will trigger an
        # exception until that rule has been run.
#        checkpoints.access_samples.get(**w)

        # Titus Brown rrightfully refers to this as magic;
        # information used to expand the pattern, is created using arbitrary Python code.
        #pattern = expand(self.pattern, sample=samples, **w)

#        rt = self2.return_read_type()

#        pattern2 = rt

#        return pattern2



import yaml
import os.path
from os import path

# checkpoint code to understand whether SE or PE workflow should be performed
class Checkpoint_GetReadType:
    def __init__(self, pattern):
        self.pattern = pattern

#    def get_samples(self):
#        # read in generated grid sample sheet
#        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
#       # get list of base_name strings
#        samples = samples_df["base_name"].tolist()
#        return samples

    def return_read_type(self):
#    if path.exists("config/read_type.yaml"):
        with open("config/read_type.yaml") as file:
            read_type_file = yaml.load(file, Loader=yaml.FullLoader)
            print(read_type_file["read_type"])
            return read_type_file["read_type"]
#    else:
#        return "PE"

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'config/grid_sample_sheet.tsv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.access_readType.get(**w)

        # Titus Brown rrightfully refers to this as magic;
        # information used to expand the pattern, is created using arbitrary Python code.
        read_type = self.return_read_type()

#        pattern = expand(self.pattern, rt=read_type, **w)
        pattern = read_type
   
        print(pattern)
        return pattern


# function checks if the read_type.yaml file already exists (created by grid_sample_preprocessing.R)
# parses yaml and then returns the read type (either SE or PE)
#def return_read_type():
#    if path.exists("config/read_type.yaml"):
#        with open("config/read_type.yaml") as file:
#            read_type_file = yaml.load(file, Loader=yaml.FullLoader)
#            print(read_type_file["read_type"])
#            return read_type_file["read_type"]
#    else:
#        return "PE"


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
#            "config/read_type.yaml"
        output:
            touch("checks/checkpoint.touch"),

    checkpoint access_readType:
        input:
#            "config/grid_sample_sheet.tsv",
            "config/read_type.yaml",
        output:
            touch("checks/checkpoint2.touch"),


    ##PE samples - GRID approach
    def get_PE_fastqs_GRID(wildcards):
        """Get raw FASTQ files based on automatically generated grid_sample_sheet.tsv."""
        rd = reads_dir
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
        s = samples_df.loc[ (wildcards.sample), ["base_name", "fq_1_end", "fq_2_end"] ].dropna()
        return [ f"{rd}/{s.base_name}{s.fq_1_end}", f"{rd}/{s.base_name}{s.fq_2_end}" ]


    ##SE samples - in the GRID approach
    def get_SE_fastqs_GRID(wildcards):
        """Get raw FASTQ files based on automatically generated grid_sample_sheet.tsv."""
        rd = reads_dir
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)
        s = samples_df.loc[ (wildcards.sample), ["base_name", "fq_1_end"] ].dropna()
        return [ f"{rd}/{s.base_name}{s.fq_1_end}" ]




#### ONLY FOR GENERAL ANALYSIS - NO STOCK MATRIX ####
elif config["approach"] == "GENERAL":
    print("GENERAL")

else:
    print("ERROR MESSAGE")
