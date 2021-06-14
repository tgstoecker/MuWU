# assign directory containing sequencing reads as specified in config.yaml
fasta = config["fasta"]
annotation = config["annotation"]
assembly_report = config["assembly_report"]


if is_gbff(annotation):

    rule GenBank_download:
        output:
            genbank="resources/annotation",
            assembly_report="resources/assembly_report.txt",
#    conda:
#        "download.yaml"
        run:
            fasta_annotation_handling(fasta, annotation, assembly_report)


    rule GenBan2gff3:




### IF USER DID NOT SUPPYLY GENBANK FILE BUT GFF(3) OR GTF  ###
if not is_gbff(annotation):
    rule fasta_annotation_download:
        output:
            genbank="resources/annotation",
            assembly_report="resources/genome.fa",
#    conda:
#        "download.yaml"
        run:
            fasta_annotation_handling(fasta, annotation, assembly_report)
