# MuWU
## Mu-Seq Workflow Utility 

- reference to publication and background information 
- special focus on python script

### Setting up the conda environment: 
Install the Python 3 version of Miniconda.
you can get it here: https://docs.conda.io/en/latest/miniconda.html

Answer yes to the question whether conda shall be put into your PATH.
For detailed options concerning conda/bioconda see:

Then, you can install Snakemake with

`conda install -c bioconda -c conda-forge snakemake`

Preparing a working directory
First, create a new directory and change into that directory in your terminal.

Download/Clone the current release of the MuWU pipeline into the directory.

With conda and the included YAML files all required software and dependencies are downloaded and prepared into conda environment during runtime of the workflow.

MuWU requires to adhere to the directory structure explained in the following.
During the workflow new directories will be created however for easy usage please copy or move your sequencing data to the RawReads directory. 
Fasta and annotation files are downloaded automatically; e.g. currently used v4 maize reference assembly and annotation:

[Zea_mays.B73_RefGen_v4.dna.toplevel.fa](ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/zea_mays/dna/)

gtf: ftp://ftp.ensemblgenomes.org/pub/plants/release-46/gtf/zea_mays/Zea_mays.B73_RefGen_v4.46.gtf.gz
gff3: ftp://ftp.ensemblgenomes.org/pub/plants/release-46/gff3/zea_mays/Zea_mays.B73_RefGen_v4.46.gff3.gz


It is also necessary to stick to the following naming scheme of the samples:
column/row; sample number; left/right
e.g.:
Col_01.1.fq; 
Col_01.2.fq; 
Row_01.1.fq; 
etc. 


### The workflow in the current release:
`snakemake --rulegraph | dot -Tsvg > rulegraph.svg`
![Alt text](./rulegraph.svg)
