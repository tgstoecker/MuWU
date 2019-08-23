# MuWU
Mu-Seq Workflow Utility 


Install the Python 3 version of Miniconda.
Answer yes to the question whether conda shall be put into your PATH.
For detailed options concerning conda/bioconda see:

Then, you can install Snakemake with

conda install -c bioconda -c conda-forge snakemake

Preparing a working directory
First, create a new directory and change into that directory in your terminal.

Download/Clone the current release of the MuWU pipeline into the directory.

The included environment.yaml file can be used to install all required software into an isolated Conda environment with a name of your choice - in the following we will call it "snakemake-MuWU":

conda env create --name snakemake-MuWU --file environment.yaml

Activating the environment
To activate the snakemake-tutorial environment, execute

conda activate snakemake-MuWU

Now you can use the installed tools and our workflow without any software dependency issues.
For detailed options of snakemake see: 

Please note that at the moment the workflow still requires to adhere to the directory structure explained in the following.
During the workflow new directories will be created however for easy usage please copy or move your sequencing data to the RawReads directory and assembly (.fa) and annotion (.gff3) to the FGS directory.

It is also recommended to stick to the following naming scheme of the samples:
column/row; sample number; left/right
e.g.:
Col_01_1.fq
Col_01_1.fq
Row_01_1.fq




In future updates I intend to increase the ease of use even more:
already planned:
- JSON file in which all options of all steps can be modified to the specific use case computional circumstance
