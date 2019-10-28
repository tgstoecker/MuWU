#####################################################################
# Identification of Mu Transposon insertion sites based on mapped   #
# Mu-Seq reads in SAM format as part of the MuSeq workflow utility  #
# MuWU                                                              #
# Tyll Stoecker & Lena Altrogge                                     #
# tyll.stoecker@uni-bonn.de                                         #
#                                                                   #
# 21.10.2019                                                        #
#                                                                   #
#####################################################################

import glob, os, re, os.path, sys, optparse, shutil
import pandas as pd
import numpy as np
import multiprocessing as mp

# Help Message

usage = "Usage: %prog --single/both -c Cores -i Inputfile -o Outputfile [Options]"

# Parse Input Parameter

parser = optparse.OptionParser(usage=usage)
parser.add_option('-i', '--input',
    dest="path", metavar="INPUT",
    help="Path to your Input Folder with mapped Mu-Seq reads in SAM Format")
parser.add_option('-o', '--output',
    metavar="OUTPUTFILE",dest="OutputFile",
    help="The name of the Output File")
#parser.add_option('--all',dest="all",
#    action='store_true', help="Output all insertions")
parser.add_option('--single', action='store_true',
    dest="single", help="Output only Single line insertions")
parser.add_option('--both', dest="both", action='store_true',
    help="Output all AND Single line insertions in two seperate files")
parser.add_option('--threshold', default=2,
    dest="int", type="int", metavar="INT",
    help="How many reads are required to support an insertion site on each side? (default=2)")
parser.add_option('-c', '--cores', default=4,
    metavar="NUM_PROCESSORS", dest="num_processors",help="How many processing cores to use (default=4)")


options, args = parser.parse_args()

path = options.path
OutputFile = options.OutputFile
#all = options.all
single = options.single
both = options.both
threshold =options.int
#cores option has to be integer not character string!
num_processors = int(options.num_processors)

# Select the Inputfiles

os.chdir(path)
files = glob.glob("*.sam")
files = sorted(files)
print("Found the following Inputfiles:", "\n", files)
print("Started to create .tmp files:")


#read and parse files in parallel
def read_and_process_parallel(files):
    for i in files:
        tmp_files = i.replace(".sam",".csv.tmp")
        print ("Started to create:",tmp_files)
        f = pd.read_csv(i, header=None, delim_whitespace=True,low_memory=False,usecols=[0,2,3,9])
        f.columns = ['Read', 'Chr' , 'Start' , 'Seq']
        f = f[f['Start'] > 0]
        df =  pd.DataFrame(columns=['Chr','Start','End','Sample','StartReads','EndReads'])

# define empty variables

        my_list = []

# Calculate End Coordinate with length of the aligned Sequence

        for n in f.Seq:
           seqlength = len(n)
           my_list.append(seqlength)
        f['Length'] = my_list
        f['End'] = (f.Start+f.Length)-1
        f['Sample'] = tmp_files

# create empty dataframe and process each Chromosome individually

        df =  pd.DataFrame(columns=['Chr','Start','End','Sample','StartReads','EndReads'])
        for x in f.Chr.drop_duplicates():
            ftable= f.loc[f['Chr'] == x]

# Find for each End Coordinate Reads that overlap with 9 bp

            for End in ftable.End.drop_duplicates():
                ftableStart = ftable.loc[ftable['Start'] == End-8]
                if len(ftableStart.index) >= threshold:
                    ftableStart = ftableStart[['Chr','Start','Sample']]
                    ftableStart['StartReads'] = len(ftableStart.index)
                    ftableStartEnd = ftable.loc[ftable['End'] == End]
                    if len(ftableStartEnd.index) >= threshold:
                        ftableStart['End'] = End
                        ftableStart['EndReads'] = len(ftableStartEnd.index)
                        ftableStart = ftableStart.drop_duplicates()
                        ftableStartfinish = ftableStart[['Chr','Start','End','Sample','StartReads','EndReads']]
                        df = df.append(ftableStartfinish)
                else:
                        continue
        df.to_csv(tmp_files,index = False)
    return;

#multiprocessing options

pool = mp.Pool(num_processors)
chunks = [files[i::num_processors] for i in range(num_processors)]
result = pool.map(read_and_process_parallel, chunks)
pool.close

#define all tmp files outside of the loop
extension = 'csv.tmp'
all_tmp_files = [i for i in glob.glob('*.{}'.format(extension))]

#combine all files in the list

print("Concatenation of all temporary files...")
combined_csv = pd.concat([pd.read_csv(f) for f in all_tmp_files])

#remove file extensions from sample column in the OutputFile (_dedup.csv.tmp gets removed in the dataframe = 14)

print("Removal of file extensions...")
combined_csv['Sample'] = [x[:-14] for x in combined_csv['Sample']]

#perform column sorting - 1.sample, 2.chromosome, 3.start coordinate

print("Sorting by sample, chromosome and start position...")
#sort_by_sample_chromosome_start = combined_csv.sort_values(['Sample','Chr', 'Start'])
combined_csv = combined_csv.sort_values(['Sample','Chr', 'Start'])
combined_csv_option_both = combined_csv.sort_values(['Sample','Chr', 'Start'])

#if the --single or --both options are chosen:

if options.single or options.both :
    print ('Identifying Single Line Insertions...')
    MuTableSingle = pd.DataFrame(columns=['Chr','Start','End','Sample','StartReads','EndReads'])
    combined_csv['combined'] = list(zip(combined_csv.Chr, combined_csv.Start, combined_csv.End))
    fRow=combined_csv[combined_csv['Sample'].str.match('Row')]
    fCol=combined_csv[combined_csv['Sample'].str.match('Col')]
    intersection=set(fRow.combined).intersection(fCol.combined)

# turn set object into list
    intersection_list = list(intersection)

    def split_list(alist, wanted_parts=1):
        length = len(alist)
        return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
                 for i in range(wanted_parts) ]

    intersection_chunks = split_list(intersection_list, wanted_parts=num_processors)
    
    def f(args):
        ''' Perform computation and write
        to separate file for each '''
        x = args[0]
        fname = args[1]
        df2 = pd.DataFrame(columns=['Chr','Start','End','Sample','StartReads','EndReads','combined'])
        df3 = pd.DataFrame(columns=['Chr','Start','End','Sample','StartReads','EndReads'])
        for i in x:
            if len(fCol.loc[fCol['combined'] == i].index) == 1 and len(fRow.loc[fRow['combined'] == i].index) ==1:
                df2Col = fCol.loc[fCol['combined'] == i]
                df2Row = fRow.loc[fRow['combined'] == i]
                df2 = df2.append(df2Col[['Chr','Start','End','Sample','StartReads','EndReads']], sort=True)
                df2 = df2.append(df2Row[['Chr','Start','End','Sample','StartReads','EndReads']], sort=True)
            else:
                continue
        df3 = df3.append(df2, sort=True)
        df3 = df3[['Chr','Start','End','Sample','StartReads','EndReads']]
        header = ['Chr','Start','End','Sample','StartReads','EndReads']
        with open(fname, 'w') as fout:
            df3.to_csv(fout, index=False, columns = header)

    if __name__ == '__main__':
        # Each sublist is a combination
        # of arguments - number and temporary output
        # file name
        x = range(len(intersection_chunks))
        names = [str(y) + '.single.tmp.csv' for y in x]
        args = list(zip(intersection_chunks,names))

        p = mp.Pool(num_processors)
        p.map(f, args)

        p.close()
        p.join()

#define all tmp files outside of the loop
    extension = 'single.tmp.csv'
    all_single_tmp_files = [i for i in glob.glob('*.{}'.format(extension))]

#combine all files in the list

    print("Concatenation of all temporary files...")
    combined_single_csv = pd.concat([pd.read_csv(i) for i in all_single_tmp_files])

#perform column sorting - 1.sample, 2.chromosome, 3.start coordinate

    print("Sorting by sample, chromosome and start position...")
    combined_single_csv = combined_single_csv.sort_values(['Sample','Chr', 'Start'])

#export to csv, create Output file and move to new directory
if not os.path.exists('../MuSeq_table'):
    os.makedirs('../MuSeq_table')


if options.single:
    combined_single_csv.to_csv('SLI-'+OutputFile, index = False)
#    shutil.move('SLI-'+OutputFile, "../MuSeq_table/",'SLI-'+OutputFile)

if options.both:
    combined_single_csv.to_csv('SLI-'+OutputFile, index = False)
    shutil.move('SLI-'+OutputFile, "../MuSeq_table_final/",'SLI-'+OutputFile)
#and
    combined_csv_option_both.to_csv(OutputFile,index = False)
    shutil.move(OutputFile, "../MuSeq_table_final/",OutputFile)

if not options.single and not options.both:
    combined_csv.to_csv(OutputFile,index = False)
#    shutil.move(OutputFile, "../MuSeq_table/",OutputFile)

#delete all tmp_files

print("Deleting temporary files...")

if 'all_tmp_files' in globals():
    for tmp_file in all_tmp_files:
        try:
            os.remove(tmp_file)
        except:
            print("Error while deleting file : ", tmp_file)

if 'all_single_tmp_files' in globals():
    for tmp_file in all_single_tmp_files:
        try:
            os.remove(tmp_file)
        except:
            print("Error while deleting file : ", tmp_file)

#Final message

if options.single:
    print("Finished creating output file :", 'SLI-'+OutputFile, "in directory ../MuSeq_table")

if options.both:
    print("Finished creating output files :", 'SLI-'+OutputFile, "and", OutputFile, "in directory ../MuSeq_table")

if not options.single and not options.both:
    print("Finished creating output file :", OutputFile, "in directory ../MuSeq_table")
