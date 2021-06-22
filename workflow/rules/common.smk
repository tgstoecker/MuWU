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
         

    ## regardless of PE or SE for fastqc we need simply all samples to be returned
    def get_fastqs_fastqc_1_GRID(wildcards):
        rd = reads_dir
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)

        s = samples_df.loc[ (wildcards.sample), ["base_name", "fq_1_end"] ].dropna()
        return [ f"{rd}/{s.base_name}{s.fq_1_end}" ]

    def get_fastqs_fastqc_2_GRID(wildcards):
        rd = reads_dir
        samples_df = pd.read_csv("config/grid_sample_sheet.tsv", dtype=str, sep="\t").set_index(["base_name"], drop=False)

        # if fq_2_end is a column we continue in paired-end mode
        if 'fq_2_end' in samples_df.columns:
            s = samples_df.loc[ (wildcards.sample), ["base_name", "fq_2_end"] ].dropna()
            return [ f"{rd}/{s.base_name}{s.fq_2_end}" ]


#### ONLY FOR GENERIC ANALYSIS - NO STOCK MATRIX ####
elif config["approach"] == "GENERIC":
#    print("GENERIC")

    samples = pd.read_csv("config/samples.tsv", dtype=str, sep="\t").set_index(["sample"], drop=False)
    SAMPLES = samples['sample'].to_list()

    def is_single_end_GENERIC_experiment(SAMPLES):
        fq2_present = pd.isnull(samples.loc[(SAMPLES), "fq2"]).to_list()
        if not any(fq2_present):
            return False 
        else:
            return True


    def is_single_end_GENERIC_sample(sample):
        """Determine whether unit is single-end."""
        fq2_present = pd.isnull(samples.loc[(sample), "fq2"])
        if isinstance(fq2_present, pd.core.series.Series):
            # if this is the case, get_fastqs cannot work properly
            raise ValueError(
                f"Multiple fq2 entries found for one sample {sample}.\n"
                "This is most likely due to a faulty samples.tsv file, e.g. "
                "Try checking your samples.tsv for duplicates."
            )
        return fq2_present


    def get_fastqs_GENERIC(wildcards):
        rd = reads_dir
        """Get raw FASTQ files from unit sheet."""
        if is_single_end_GENERIC_sample(wildcards.sample):
            s = samples.loc[ (wildcards.sample), ["fq1"] ].dropna()
            return [ f"{rd}/{s.fq1}" ]
        else:
            u = samples.loc[ (wildcards.sample), ["fq1", "fq2"] ].dropna()
            return [ f"{rd}/{u.fq1}", f"{rd}/{u.fq2}" ]


    def get_fastqs_fastqc_1_GENERIC(wildcards):
        rd = reads_dir
        """Get raw FASTQ files from sample sheet."""
        s = samples.loc[ (wildcards.sample), ["fq1"] ].dropna()
        return [ f"{rd}/{s.fq1}" ]

    def get_fastqs_fastqc_2_GENERIC(wildcards):
        rd = reads_dir
        """Get raw FASTQ files from sample sheet."""
        s = samples.loc[ (wildcards.sample), ["fq2"] ].dropna()
        return [ f"{rd}/{s.fq2}" ]


else:
    print("ERROR MESSAGE")
