import os
import posixpath

from urllib.error import URLError
from urllib.request import urlopen
from urllib.request import retrieve


#pseudocode ideas
if annotation ending = .gbff skip anything fasta related


def is_url(file):
    return (
        file.startswith("http")
        or file.startswith("ftp")
        or file.startswith("sftp")
    )

def get_file(file):
    if not is_url(file):
        return file
    else:
        

def find_extension(file, extensions=[".gff", ".gff3", ".gbff", ".gtf"]):
    for ext in extensions:
        if file.endswith("{}".format(ext)):
            return file
        else:
            try:
                urlopen(file + script)
                return file + script
            except URLError:
                continue




#if config option download fasta:
rule fasta_annotation_download:
    output:
        fa="FGS/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",
#        gtf="FGS/Zea_mays.B73_RefGen_v4.50.gtf",
#        gff3="FGS/Zea_mays.B73_RefGen_v4.50.gff3",
    params: #under these locations in the config file the links to the respective files are listed
        fasta=config["fasta"],
        gtf=config["gtf"],
        gff3=config["gff3"],
    conda:
        "download.yaml"
    shell:
        """
        wget -O - {params.fasta} | gunzip -c > {output.fa}
#        wget -O - {params.gtf} | gunzip -c > {output.gtf}
#        wget -O - {params.gff3} | gunzip -c > {output.gff3}
        """


#if config option download annotation:
rule fasta_annotation_download:
    output:
#        fa="FGS/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",
#        gtf="FGS/Zea_mays.B73_RefGen_v4.50.gtf",
        gff3="FGS/Zea_mays.B73_RefGen_v4.50.gff3",
    params: #under these locations in the config file the links to the respective files are listed
        fasta=config["fasta"],
        gtf=config["gtf"],
        gff3=config["gff3"],
    conda:
        "download.yaml"
    shell:
        """
#        wget -O - {params.fasta} | gunzip -c > {output.fa}
#        wget -O - {params.gtf} | gunzip -c > {output.gtf}
        wget -O - {params.gff3} | gunzip -c > {output.gff3}
        """

#works on gff or gtf
rule extract_gene_transcript_tables:
#gene
gffread Zea_mays.B73_RefGen_v4.50.gff3 --keep-genes --table @geneid,@chr,@start,@end -o test
#transcript # perhaps add final blacklist to remove prefixes "transcript:"
#gffread Zea_mays.B73_RefGen_v4.50.gff3 --table @id,@chr,@start,@end -o test

#gffread version 0.12.1



#if gbff format - must provide assembly_report as well

#if link is provided (either case) download file
#alternatively use file as given path onf file system


#mamba install -c anaconda biopython=1.78

#MIT license biocode code repurposed and altered here / under
python convert_genbank_to_gff3.py --with_fasta 
-i GCA_000005005.6_B73_RefGen_v4_genomic.gbff 
-o GCA_000005005.6_B73_RefGen_v4_genomic.gff3

#fasta is appended and we can split on "##FASTA"

#csplit 8.30
#csplit -s -z GCA_000005005.6_B73_RefGen_v4_genomic.gff3 /##FASTA/ '{*}'
#xx00 = gff3
#xx01 = FASTA

mamba install -c conda-forge coreutils=8.31

#suppress-matched will remove the ##FASTA line on which we split
csplit -s --suppress-matched -z ../GCA_000005005.6_B73_RefGen_v4_genomic.gff3 /##FASTA/ '{*}' && 
mv xx00 annotation.gff3 && 
mv xx01 genome.fasta && 
mv annotation.gff3 genome.fasta ../test_dir2/

#use report.txt to rename chromosome column
grep -v '#' GCA_000005005.6_B73_RefGen_v4_assembly_report.txt | 
awk -F "\t" '{print $5,$1}' OFS="\t" | 
awk -F "\t" '{sub(/\..*$/,"",$1)}1' OFS="\t" > file_used_for_renaming


#using seqkit to rename fasta headers
mamba install -c bioconda seqkit=0.16.1 
seqkit replace -p "(.+)" -r '{kv}' -k ../file_used_for_renaming headers_test | less

#with awk and sed exchange chromosome names in gff3 file
awk -F "\t" 'NR==FNR{t[$1]=$2} NR!=FNR && FNR>1{$1=t[$1];print}' OFS="\t" ../file_used_for_renaming annotation.gff3 | 
sed '1 i\##gff-version 3' > changed_annotation


#report file may contain NAs if genome doesnt contain Mt information etc.
#that's why we need sth- like a left join to only try:
#get name from genome/annoation; check if it is inside report; rename it to first column

#we get 
#bioperl approach with bp_genbank2gff3.pl - sadly much too slow 
#mamba install -c biobuilds perl=5.22
#mamba install -c bioconda perl-bioperl
#mamba install -c bioconda perl-cgi


