# MuWU
## Mu-Seq Workflow Utility 

- Automated workflow for the identification of Mutator insertion sites used in the creation of the BonnMu resource

### Setting up the conda environment: 
Install the Python 3 version of Miniconda.
you can get it here: https://docs.conda.io/en/latest/miniconda.html

Answer yes to the question whether conda shall be put into your PATH environment variable.
For detailed options concerning conda/bioconda see:

Then, you can install Snakemake with

`conda install -c bioconda -c conda-forge snakemake`

### Preparing a working directory:

Download/Clone the current release of the MuWU pipeline.

With conda and the included YAML files all required software and dependencies are downloaded and prepared into conda environment during runtime of the workflow.

MuWU requires to adhere to the directory structure explained in the following.
During the workflow new directories will be created, however for easy usage please copy or move your sequencing data to the RawReads directory. 
Fasta and annotation files are downloaded automatically from ensembl; e.g. currently used: v4 maize reference assembly and annotations:

- Zea_mays.B73_RefGen_v4.dna.toplevel.fa

- Zea_mays.B73_RefGen_v4.46.gtf

- Zea_mays.B73_RefGen_v4.46.gff3


It is also necessary to stick to the following naming scheme of the samples:
column/row; sample number; left/right
e.g.:
Col_01.1.fq; 
Col_01.2.fq; 
Row_01.1.fq; 
etc. 

### Starting MuWU:
Change thread options for individual rules in the config.yaml file. Then specifiy overall threads and start MuWU via:
`snakemake --cores xx --use-conda



### The workflow in the current release:
`snakemake --rulegraph | dot -Tsvg > rulegraph.svg`
![Alt text](./rulegraph.svg)
