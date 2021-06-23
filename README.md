# MuWU - Mu-Seq Workflow Utility [![Snakemake](https://img.shields.io/badge/snakemake-=6.4.1-brightgreen.svg)](https://snakemake.readthedocs.io)

<br>

<img align="right" src="TSD_mapping.png" width=300 >


- Automated workflow for the identification and annotation of transposable element insertion sites originally developed for the BonnMu resource and *Mutator* transposons in particular 
- MuWU is able to **detect any kind of TE insertion event** as long as target site duplications (TSDs) are created by its integration and the TSD length is known

<br>  

# :control_knobs: Two modes - GRID and GENERIC:

**GRID**
- Requires as input reads in grid design as outlined e.g. by McCarty et al. 2013; Liu et al. 2016; Marcon et al. 2020 
- Differentiates between heritable germinal insertions and somatic insertions and annotates both sets

**GENERIC**
- Does not require sequencing reads of special experimental design
- Identifies & annotates all insertions of the particular TE


<br>  

# :gear: Options

config.yaml, bla

<br>  

# :arrow_double_down: Download & Setup

There are 2 ways of using MuWU:  
* via cloning this repo and then using conda installation of necessary software at runtime
* a singularity container, which includes a and requires no further downloads except for the container itself  


## Option 1. Cloning of this repo and download/installation of software at runtime
### Step 1 - Set up conda and snakemake: 
Install the Python 3 version of Miniconda.
you can get it here: https://docs.conda.io/en/latest/miniconda.html

Answer yes to the question whether conda shall be put into your PATH environment variable.

Then, you can install mamba and Snakemake with

`conda install -c conda-forge -c bioconda mamba snakemake=6.4.1`  

### Step 2 - Preparing the working directory:

Clone the current release of the MuWU pipeline.

`git clone https://github.com/tgstoecker/MuWU.git`

With conda and the included YAML files all required software and dependencies are downloaded and prepared into conda environment during runtime of the workflow.


## Option 2. Singularity container
### Step 1 - Set up Singularity on your system: 
Install the Python 3 version of Miniconda.
you can get it here: https://docs.conda.io/en/latest/miniconda.html

Answer yes to the question whether conda shall be put into your PATH environment variable.

Install mamba:
`conda install -c conda-forge mamba`

Then, you can install Singularity (3.6.1) with
`mamba install -c conda-forge singularity=3.6.1`  
  
Alternatively install Singularity based on these instructions: https://singularity.lbl.gov/install-linux  
  
### Step 2 - set up the container to run the workflow  ####

**Download the MuWU-1.1.sif file, hosted here:**  https://uni-bonn.sciebo.de/s/LsmuNDuEeA0sUer  
  
Create a sandbox from the .sif file:  
This can take a while (on Intel(R) Xeon(R) CPU E5-2690 v4@ 2.60GHz roughly 30 min!), since the sandbox will be over 20Gb in size.  
It is might be necessary to set SINGULARITY_TMPDIR to a particular (or newly created) tmp directroy as singularity on some systems uses `/tmp` directory as standard while building. This can lead to storage errors if the space is limited by your sysadmin.  
Easy workaround - set SINGULARITY_TMPDIR to a directory where space is plenty:  
`export SINGULARITY_TMPDIR=/path/to/where/tmp/should/be`  
  
`singularity build --sandbox MuWU-example MuWU-example.sif`  
  
Access the sandbox:  
`singularity shell --writable --no-home MuWU-example`
  
Once "inside", navigate to the MuWU directory  
`cd /MuWU/`  
  
Activate conda environment (snakemake is already installed):  
`source activate snakemake`  

<br>  

# :beginner: Usage & Output

## Required input files

As described under options, control of parameters and inputs is inside `config/config.yaml` - for more details for all options please see that file.

Both the GRID & GENERIC methods require:
1. Reference sequence & annotation for the species in question 
   - MuWU can handle both file paths as well URL links (will download files in the later case automatically)
   - Files can be either unpacked or gzipped
   - We currently support gff3, gtf and genbank (.gbff & .dat) as annotation formats
     - In case of a GenBank annotation we also demand the corresponding "assembly_report.txt" to be supplied in order to cirrectly rename the chromosomes/scaffolds
2. SE or PE sequencing reads (unpacked or gzipped) (best: enriched for insertions and with a primer/adapter PCR approach yielding starts/ends with TSD sequence after trimming)
3. File describing samples (**this differs between the methods!**)
   - GRID: an excel table needs to be supplied under `config/stock_matrix/` (example provided)
   - GENERIC: appropiately modified `config/samples.tsv` file
   - -> In both cases the file is used to infer the complete structure of the workflow and SE/PE type of the reads
   - -> In the GRID method it is important that the base name of the fastq file/s match its corresponding entry in the the stock matrix table


## Once everything is set up - run the workflow:  

Check the workflow (dryrun; testbuild of DAG):  
`snakemake --use-conda --cores 24 --conda-prefix conda_envs -np`
  
Run the workflow:  
`snakemake --use-conda --cores 24 --conda-prefix conda_envs`  

`--conda-prefix conda_envs` will look for/install the environments in a directory called `conda_envs/`.  
This is especially important if you should use the singularity container. Here the main software and test folder all have their respective environments installed in such a directory. If you omit this parameter snakemake/conda will try to download & install all required software which the container was specifically build for to circumvent.


## Output  
Besides a Final outputs are generated in the directory 
The main output files are:  

1. MultiQC HTML output (open in browser):  
```
/MuWU/multiqc/multiqc.html
```  
  
2. (Annotated) Insertion tables under `MuWU/results/insertions_table_final/`  
```
all_identified_insertions_annotated.csv

#! only 
l_identified_insertions_annotated.csv

```

<br>  


### Starting MuWU:
Change thread options for individual rules in the config.yaml file.  
Then specifiy overall threads and start MuWU via:  
`snakemake --cores xx --use-conda`  
<br>  

# :heavy_check_mark: Tests
We have included several tests to demonstrate MuWU' functionality and broad range of usability.

- 1
- 2
- 3

<br>  

# :framed_picture: Visualized
### GRID workflow (as snakemake rulegraph):
`snakemake --rulegraph | dot -Tsvg > rulegraph.svg`

<br>  

![Alt text](./rulegraph_GRID.svg)
