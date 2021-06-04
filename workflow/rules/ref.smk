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


rule extract_gene_transcript_tables:
#gene
gffread Zea_mays.B73_RefGen_v4.50.gff3 --keep-genes --table @geneid,@chr,@start,@end -o test
#transcript # perhaps add final blacklist to remove prefixes "transcript:"
gffread Zea_mays.B73_RefGen_v4.50.gff3 --table @id,@chr,@start,@end -o test

#gffread version 0.12.1
