import os
import posixpath
import gzip
import shutil

from urllib.request import urlretrieve
from urllib.request import urlcleanup


extensions=[".gff", ".gff3", ".gbff", ".gtf", ".dat", ".gff.gz", ".gff3.gz", ".gbff.gz", ".gtf.gz", ".dat.gz"]
genbank=[".gbff", ".dat", ".gbff.gz", ".dat.gz"]


def is_url(file):
    return (
        file.startswith("http")
        or file.startswith("ftp")
        or file.startswith("sftp")
    )

def is_valid_annotation_file(annotation):
    for ext in extensions:
        if annotation.endswith("{}".format(ext)):
            return True

def get_file_ext(annotation):
    for ext in extensions:
        if annotation.endswith("{}".format(ext)):
            return 

def is_gbff(file):
    for ext in genbank:
        if file.endswith("{}".format(ext)):
            return True

def is_gzipped(file):
    return(
        file.endswith(".gz")
    )

def handle_fasta(fasta, annotation):
    if not is_gbff(annotation):
        if is_url(fasta):
            print("Downloading fasta")
            if is_gzipped(fasta):
                urlretrieve(fasta,'resources/genome.fa.gz')
                print("Unpacking fasta")
                with gzip.open('resources/genome.fa.gz', 'rb') as f_in:
                    with open('resources/genome.fa', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            else:
                urlretrieve(fasta,'resources/genome.fa')
        elif not is_url(fasta):
            print("Gunzipping or linking fasta file already on file system to resources/ dir")
            if is_gzipped(fasta):
                with gzip.open(fasta, 'rb') as f_in:
                    with open('resources/genome.fa', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            else:
                os.symlink(os.getcwd() + "/" + fasta, "resources/genome.fa")

def gunzip_annotation(annotation_file):
    with gzip.open(annotation_file, 'rb') as f_in:
        with open(annotation_file.replace(".gz", ""), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def handle_annotation(annotation):
    urlcleanup()
    if is_valid_annotation_file(annotation) and not is_gzipped(annotation):
        if is_url(annotation):
            print("Downloading annotation")
            urlretrieve(annotation, 'resources/annotation')
        elif not is_url(annotation):
            printing("Linking annotation")
            os.symlink(annotation, 'resources/annotation' + get_file_ext(annotation))
    if is_valid_annotation_file(annotation) and is_gzipped(annotation):
        if is_url(annotation):
            print("Downloading annotation")
            urlretrieve(annotation, 'resources/annotation.gz')
        elif not is_url(annotation):
            print("Linking annotation")
            os.symlink(annotation, "resources/annotation.gz")
        print("Unpacking annotation")
        gunzip_annotation('resources/annotation.gz')
        print("Done!")

def fasta_annotation_handling(fasta, annotation, assembly_report):
    # check if annotation is valid file - based on file extensions
    if is_valid_annotation_file(annotation):
        print("Annotation file type is supported - continuing:")
        # dealing with the fasta file
        # download/link fasta to resources/ dir
        # we only download/link/use fasta if annotation is NOT in GenBank format
        handle_fasta(fasta, annotation)
        # dealing with the annotation file
        handle_annotation(annotation)
        if is_gbff(annotation):
            print("Supplied annotation in GenBank format - also getting assembly report.")
            if is_url(assembly_report):
                urlretrieve(assembly_report, 'resources/assembly_report.txt')
            else:
                os.symlink(assembly_report, "resources/assembly_report.txt")
        print("All done!")
    else:
        print("Unrecognized annotation file type!")
        print("Valid files are:")
        print(extensions)


#fasta_annotation_handling(fasta, annotation, report)

