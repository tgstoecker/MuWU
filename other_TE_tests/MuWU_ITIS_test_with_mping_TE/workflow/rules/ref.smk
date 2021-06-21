# assign directory containing sequencing reads as specified in config.yaml
fasta = config["fasta"]
annotation = config["annotation"]
assembly_report = config["assembly_report"]


if is_gbff(annotation):

    rule GenBank_path:
        output:
            # we rename the standardized "annotation" output here in order to not run into issues with inprecise rules later
            genbank="resources/annotation_gb",
            assembly_report="resources/assembly_report.txt",
#        conda:
#            "../envs/download.yaml"
        run:
            fasta_annotation_handling(fasta, annotation, assembly_report)
            os.rename(r'resources/annotation',r'resources/annotation_gb')


    rule GenBank2gff3:
        input:
            "resources/annotation_gb"
        output:
            "resources/annotation.gff3"
        log:
            "logs/GenBank2gff3/GenBank2gff3/convert_genbank_to_gff3.log"
        conda:
            "../envs/genbank_to_gff3.yaml"
        #MIT license biocode code repurposed and altered here / under scripts/
        shell:
            "python workflow/scripts/convert_genbank_to_gff3.py --with_fasta -i {input} -o {output} 1>{log} 2>>{log}"


    #fasta is appended and we can split on "##FASTA"
    rule split_annotation_fasta:
        input:
            "resources/annotation.gff3"
        output:
            split_annotation = "resources/split_annotation.gff3", 
            fa = "resources/split_genome.fa",
        conda:
            "../envs/coreutils.yaml"
        # {{*}} = results inside shell to {*}
        #suppress-matched will remove the ##FASTA line on which we split
        shell:
            "csplit -s --suppress-matched -z {input} /##FASTA/ '{{*}}' && "
            "mv xx00 {output.split_annotation} && "
            "mv xx01 {output.fa} "


    rule create_renaming_scheme:
        input:
            "resources/assembly_report.txt"
        output:
            "resources/renaming_scheme.txt"
        conda:
            "../envs/renaming.yaml"
        # note the double backslash in front of t - "\\t" becomes "\t" in shell
        # we create a table with a genbank accession - sequence name AND a refseq accession - sequence name relationship
        # both relationships are after one another - long format; this works with both fasta as well as annotation renaming and is thus flexible
        # to whatever one of the two options the user has decided to get the GenBank file from 
        shell:
            """
            grep -v '#' {input} |
            awk -F "\\t" '{{print $7,$1}}' OFS="\t" |
            awk -F "\\t" '{{sub(/\..*$/,"",$1)}}1' OFS="\t" > {output} &&
            grep -v '#' {input} |
            awk -F "\\t" '{{print $7,$1}}' OFS="\t" |
            awk -F "\\t" '{{sub(/\..*$/,"",$1)}}1' OFS="\t" >> {output}
            """

    rule rename_fasta_headers:
        input:
            scheme = "resources/renaming_scheme.txt",
            fa = "resources/split_genome.fa",
        output:
            "resources/genome.fa"
        conda: "../envs/seqkit.yaml"
        shell:
            "seqkit replace -p '(.+)' -r '{{kv}}' -k {input.scheme} {input.fa} > {output}"


    rule rename_annotation_chromosomes:
        input:
            scheme = "resources/renaming_scheme.txt",
            split_annotation = "resources/split_annotation.gff3",
        output:
            "resources/annotation"
        conda: "../envs/coreutils.yaml"    
        # important note see: rule create_renaming_scheme
        shell:
            """
            awk -F '\\t' 'NR==FNR{{t[$1]=$2}} NR!=FNR && FNR>1{{$1=t[$1];print}}' OFS='\\t' {input.scheme} {input.split_annotation} |
            sed '1 i\\##gff-version 3' > {output}
            """


### IF USER DID NOT SUPPYLY GENBANK FILE BUT GFF(3) OR GTF  ###
if not is_gbff(annotation):
    rule fasta_gff_gtf_path:
        output:
            fa = "resources/genome.fa",
            gff_gtf = "resources/annotation",
#    conda:
#        "download.yaml"
        run:
            fasta_annotation_handling(fasta, annotation, assembly_report)


### regardless of input files - create a common final table format and a bowtie2 index for the reference genome ###
rule create_annotation_tables:
    input:
        "resources/annotation"
    output:
        "resources/final_annotation_table"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --keep-genes --table @geneid,@chr,@start,@end -o {output}"


rule bowtie2_index:
    input:
        reference = "resources/genome.fa"
    output:
        multiext(
            "resources/genome",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log"
    params:
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "file:workflow/builds/MuWU_bowtie2/bowtie2_build"
    




