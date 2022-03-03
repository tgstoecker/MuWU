#for now this rule expects gzipped, paired data and GRID design
#I can comeback sometime and implement it dynamically for every scenario

rule read_te_typing:
    input:
        "rawreads/{sample}_{paired}.fq.gz",
#        check="checks/read_type.check",
#        fastq_files=get_fastqs_GRID,
    output:
        "results/te_typing/pre_sorting/{sample}/{paired}/{te_type}.txt",
    params:
    #noice, access to nested yaml via lambda function
    # everything stays tidy
        type_seq=lambda w: config["TE_types"][w.te_type],
        type=lambda w: config["TE_types"],
#    log:
#       "logs/cutadapt/{sample}.log"
    threads: config["threads_te_typing_seqkit"]
    conda: "../envs/te_typing.yaml"
    shell:
#        "echo {params.type_seq}"
        "zcat {input} | "
        "seqkit grep --by-seq -m 0 -j {threads} -p {params.type_seq} | " 
        "seqkit seq --name | " 
        "sed 's/\s.*$//' | "
        """awk '{{print $1,"{wildcards.paired}","{wildcards.te_type}"}}' OFS="\\t" > {output}"""


rule merging_read_te_typing:
    input:
        expand("results/te_typing/pre_sorting/{sample}/{paired}/{te_type}.txt", sample=SAMPLES, paired=PAIRED, te_type=TE_TYPES)
#    params:
#        type_seq=lambda w: config["TE_types"][w.te_type],
#        type=lambda w: config["TE_types"],
    output:
        "results/te_typing/pre_sorting/{sample}/{sample}_te_types_merged.tsv"
    shell:
         """
         cat results/te_typing/pre_sorting/{wildcards.sample}/1/* results/te_typing/pre_sorting/{wildcards.sample}/2/* > {output} &&
         sed  -i '1i Name\\tStrand\\tType' {output}
         """


####### Annotation with te type/s  #######

rule te_typing_annotation:
    input:
#        germinal_annotated="results/insertions_table_final/germinal_identified_insertions_annotated.csv",
#        germinal="results/insertions_table_final/germinal_identified_insertions.csv",
#        all_annotated="results/insertions_table_final/all_identified_insertions_annotated.csv",
#        all="results/insertions_table_final/all_identified_insertions.csv",
        insertion_table="results/insertions_table_final/{insertion_table}.csv",
    output:
        "results/insertions_table_final_te_typed/headers_strand_1_uncategorized_{insertion_table}.csv",
        "results/insertions_table_final_te_typed/headers_strand_2_uncategorized_{insertion_table}.csv",
    params:
        samples=lambda wildcards: list(config["SAMPLES"]),
        all_types=lambda wildcards: list(config["TE_types"].keys()), 
        insertion_table_name=lambda wildcards: wildcards.insertion_table,
    conda: "../envs/annotation.yaml"
    script:
        "../scripts/te_type_annotation.R"
#    shell:
#        """
#        {params.samples}
#        {params.all_types}
#        """


###### once done with annotation  ###########
#extract all reads from fqs that correspond to uncategorized insertions

rule get_uncategorized_ins_reads_1:
    input:
        reads="rawreads/{sample}_1.fq.gz",
        unc_headers="results/insertions_table_final_te_typed/headers_strand_1_uncategorized_{insertion_table}.csv",
    output:
        "results/te_typing/uncategorized/{insertion_table}/{sample}/1/unc.fa",
    threads: config["threads_te_typing_seqkit"]
    conda: "../envs/te_typing.yaml"
    shell:
        """
        zcat {input.reads} |
        seqkit grep -j {threads} --by-name -r -f {input.unc_headers} | 
        seqkit grep -j {threads} --by-seq -m 4 -p ATAATGGCAATTATCTC |
        seqkit seq --seq-type DNA --reverse --complement |
        seqkit fq2fa > {output}
        """

rule get_uncategorized_ins_reads_2:
    input:
        reads="rawreads/{sample}_2.fq.gz",
        unc_headers="results/insertions_table_final_te_typed/headers_strand_2_uncategorized_{insertion_table}.csv",
    output:
        "results/te_typing/uncategorized/{insertion_table}/{sample}/2/unc.fa",
    threads: config["threads_te_typing_seqkit"]
    conda: "../envs/te_typing.yaml"
    shell:
        """
        zcat {input.reads} |
        seqkit grep -j {threads} --by-name -r -f {input.unc_headers} |
        seqkit grep -j {threads} --by-seq -m 4 -p ATAATGGCAATTATCTC |
        seqkit fq2fa > {output}
        """

rule merge_uncategorized_ins_reads:
    input:
#        expand("results/te_typing/uncategorized/{sample}/{{paired}}/unc.fa", sample=SAMPLES),
        one=expand("results/te_typing/uncategorized/{sample}/1/unc.fa", sample=SAMPLES),
        two=expand("results/te_typing/uncategorized/{sample}/2/unc.fa", sample=SAMPLES),
    output:
#        "results/te_typing/uncategorized/merged/{paired}/merged_{paired}_unc.fa",
        one="results/te_typing/uncategorized/merged/1/merged_1_unc.fa",
        two="results/te_typing/uncategorized/merged/2/merged_2_unc.fa",
    conda: "../envs/te_typing.yaml"
    shell:
        """
        cat {input.one} > {output.one}
        cat {input.two} > {output.two}
        """


#create faidx index
rule index_categorized_ins_reads:
    input:
        "results/te_typing/uncategorized/merged/{paired}/merged_{paired}_unc.fa",
    output:
        "results/te_typing/uncategorized/merged/{paired}/merged_{paired}_unc.fa.fai"
#    params:
#        p=lambda w: config["PAIRED"]
    conda: "../envs/te_typing.yaml"
    shell:
        """
        samtools faidx {input}
        """

#create bed file with information where in the read the more or less conserved motif is situated
#use the same mismatch value (for now = 4) & mostly conserved motif ATAATGGCAATTATCTC
rule locate_motif:
    input:
        fa="results/te_typing/uncategorized/merged/{paired}/merged_{paired}_unc.fa",
        fai="results/te_typing/uncategorized/merged/{paired}/merged_{paired}_unc.fa.fai"
    output: 
        "results/te_typing/uncategorized/merged/{paired}/merged_{paired}_unc_motif_info.bed"
    conda: "../envs/te_typing.yaml"
    shell:
        """
        cat {input.fa} | seqkit locate -m 4 -p ATAATGGCAATTATCTC > {output}
        """


###!depending on motif choice!
#modify the bed file so that the start coordinate is extended to the "left" by max 12 bp
#if start coordiante is <= 12 the nset to 0
#output format is for ondex lookup: read_name:start-end
rule extend_motif:
    input:
        "results/te_typing/uncategorized/merged/{paired}/merged_{paired}_unc_motif_info.bed",
    output:
        "results/te_typing/uncategorized/merged/{paired}/{paired}_info_cut.txt",
    conda: "../envs/te_typing.yaml"
    shell:
        """awk -F "\\t" '$5 < 12 && NR>1 {{print $1":"1"-"$6; next}} NR>1 {{print $1":"$5-12"-"$6}}' {input} > {output}"""



#loop over fasta and extract subregions
#-> this would be a targeted single sequence
#seqkit faidx strand_1_final.fa E00591:404:HFCMWCCX2:7:1101:14823:10081:1-20
###conda install -c conda-forge findutils
#since we use >>; do a sanity check if an old version was already created
rule extract_unc_ins_motif_subregions:
    input:
        fa="results/te_typing/uncategorized/merged/{paired}/merged_{paired}_unc.fa",
        ext_motifs="results/te_typing/uncategorized/merged/{paired}/{paired}_info_cut.txt",
    output:
        "results/te_typing/uncategorized/merged/{paired}/{paired}_motifs_pos_header.fa",
    conda: "../envs/te_typing.yaml"
    shell:
        """
        [ -f {output} ] && rm {output}
        xargs -n 1 -I{{}} samtools faidx {input.fa} {{}} < {input.ext_motifs} >> {output}
        """


#rename header to base read + /1 or /2
#> add this point each header ends in :#-# (range in original read where the motif is)
#this is important for edge cases in which both mates of paired reads appear in the next step;
#might lead to problems, so let's be clean here
rule rename_unc_ins_motif_subregions:
    input:
        "results/te_typing/uncategorized/merged/{paired}/{paired}_motifs_pos_header.fa",
    output:
        "results/te_typing/uncategorized/merged/{paired}/renamed_{paired}_motifs_pos_header.fa",
    conda: "../envs/te_typing.yaml"
    #sed 's/\(.*\):.*/\1\\1/' {input} > {output}
    shell:
        """
        sed 's/\(.*\):.*/\\1\\\\{wildcards.paired}/' {input} > {output}
        """

### Merge/concat of _1 and _2 files here
rule concat_motif_subregions_files:
    input:
        expand("results/te_typing/uncategorized/merged/{paired}/renamed_{paired}_motifs_pos_header.fa", paired=PAIRED),
    output:
        "results/te_typing/uncategorized/merged/merged_final_motifs_pos_header.fa",
    conda: "../envs/te_typing.yaml"
    shell:
        """
        cat {input} > {output}
        """


#cluster reads - -c 1; shorter seqs if identical in their length will be "swallowed" by larger seqs
rule clustering_unc_reads_final_set:
    input:
        "results/te_typing/uncategorized/merged/merged_final_motifs_pos_header.fa",
    output:
        txt="results/te_typing/uncategorized_clustered/unc_reads_final_set.txt",
        clstr="results/te_typing/uncategorized_clustered/unc_reads_final_set.txt.clstr",
    conda: "../envs/te_typing.yaml"
    shell:
        """
        cd-hit -i {input} -o {output.txt} -c 1 -g 1 -sc 1 -sf 1 -d 80
        """ 

#cd-hit comes with reformat scripts for the output in perl
#clstr_rep column (column #5) == 1 then corresponds to representative sequence (found in ...out.txt)
#-> for this line we can look in column (#1) id to find the name of this fa entry/read
#clstr_size (column #3) == # corresponds to size of cluster
#reducing to representative seqs entries
rule reformat_clustering_results_1:
    input:
        "results/te_typing/uncategorized_clustered/unc_reads_final_set.txt.clstr",
    output:
        "results/te_typing/uncategorized_clustered/ref_unc_reads_final_set.txt.clstr",
    conda: "../envs/te_typing.yaml"
    shell:  
        """
        clstr2txt.pl {input} > {output}
        """

#re-format normal out.txt to columns
#simple solution since seqs are <=29 bp so we don't have multi-line occurences
rule reformat_clustering_results_2:
    input:
        "results/te_typing/uncategorized_clustered/unc_reads_final_set.txt",
    output:
        pre="results/te_typing/uncategorized_clustered/pre_ref_reads_final_set.txt",
        final="results/te_typing/uncategorized_clustered/final_ref_reads_final_set.txt",
    conda: "../envs/te_typing.yaml"
    shell:
        """
        cat {input} | paste - - > {output.pre}
        sed 's/^>//' {output.pre} > {output.final}
        """

#add cluster size information as third column to reformat_test_hit_out.txt
##-> create subset of reformat_test_hit_out.txt.clstr with only read and cluster size info for representative seqs/reads
## reference
#https://stackoverflow.com/questions/45167499/how-to-merge-2-tables-with-awk
rule final_cluster_sizes:
    input:
        clstr="results/te_typing/uncategorized_clustered/ref_unc_reads_final_set.txt.clstr",
        txt="results/te_typing/uncategorized_clustered/final_ref_reads_final_set.txt",
    output:
        clstr="results/te_typing/uncategorized_clustered/rep_ref_unc_reads_final_set.txt.clstr",
        final="results/te_typing/uncategorized_clustered/final_clstr_file.tsv",
    conda: "../envs/te_typing.yaml"
    shell:
        """awk -F "\\t" '$5 > 0 && NR>1 {{print $1,$3; next}} {{}}' OFS="\\t" {input.clstr} > {output.clstr} && """
        """awk 'BEGIN {{FS=OFS="\\t"}} NR==FNR {{h[$1] = $2; next}} {{print $0,h[$1]}}' {input.txt} {output.clstr} > {output.final}"""

#awk -F "\t" '$5 > 0 && NR>1{print $1,$3; next} {}' OFS="\t" reformat_test_hit_out.txt.clstr > rep_test_hit_out.txt.clstr
#awk 'BEGIN {FS=OFS="\t"} NR==FNR {h[$1] = $2; next} {print $0,h[$1]}' final_reformat_test_hit_out.txt rep_test_hit_out.txt.clstr > final_rep_clstr.txt
