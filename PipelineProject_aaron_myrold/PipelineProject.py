#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:27:05 2023

@author: aaronmyrold
"""

import os
import subprocess
from Bio import Entrez
from Bio import SeqIO
# from Bio.Blast import NCBIWWW

# determine the path of the directory this file is located in
# idea taken from here: https://www.pythonanywhere.com/forums/topic/13464/
my_env = os.path.join(os.path.dirname(__file__)) 
# my_env = '/Users/aaronmyrold/Desktop/PipelineProject/PipelineProject_aaron_myrold'
# set the current working directory to that folder so that remaining paths can function properly
os.chdir(my_env)


# PART 0 ----
# Make directories and paths
folder_names = ('1_data_raw', '2_data_clean', '3_output', '4_index', '5_data_test', '6_blast')
p_data_raw = folder_names[0]
p_data_clean = folder_names[1]
p_out = folder_names[2]
p_index = folder_names[3]
p_test = folder_names[4]
p_blast = folder_names[5]

# this is primarily due to github not saving empty directories
# need to be able to create folder structure from scratch if needed  
for i in folder_names:
    if not os.path.exists(i):
        os.makedirs(i)
        
# Make accession list
# these are used to dynamically create test files below
accessions = ['SRR5660030','SRR5660033', 'SRR5660044', 'SRR5660045']


# PART 1 ----
# Download the raw/test data
# Determine a test run or full run to determine which data set to download
raw_or_test = None                                                              
acc_val = ('raw','test')  
# while loop idea from https://stackoverflow.com/questions/37826322/get-a-specific-input-in-python                                                      
while raw_or_test not in acc_val:
    raw_or_test = input('Please type \"raw\" or \"test\" to specify a full or test run.\n')

# If we want raw data, download the data
if raw_or_test == 'raw':
    # set working directory to data_raw
    os.chdir(p_data_raw)
    # this could possibly be done with one line, but did not want to run multiple download tests
    # if its possible, would work like this:
    # os.system('fasterq-dump --threads 2 --progress SRR5660030 SRR5660033 SRR5660044 SRR5660045')
    for i in accessions:
        if not os.path.isfile(f'{i}_1.fastq.gz'):
            os.system(f'fastq-dump --gzip --split-3 --aligned {i}')
    os.chdir('..')
    
# If we want test data, create it
if raw_or_test == 'test':
    os.chdir(p_test)
    # if the test files do not exist, create them
    for i in accessions:
        if not os.path.isfile(f'{i}_1.fastq.gz'):
            os.system(f'fastq-dump -X 10000 --gzip --split-3 --aligned {i}')
    # set the p_data_raw path to test data folder
    p_data_raw = p_test
    os.chdir('..') #return to working directory


# PART 2 ----
# Here I download the reference genome (HCMV) in order to build the index for bowtie2
# I use the Entrez.efetch() method to download the viral genome
# Then I run a bowtie2 command to build the index inside of the 4_index folder
os.chdir(p_index)
Entrez.email = "amyrold@luc.edu"
# code adapted from https://www.biostars.org/p/261774/
# code originally from biopython cookbook, wasn't sure which to link so I linked the site I found it from.
filename = "NC_006273.2.fasta"
if not os.path.isfile(filename): # To prevent overwriting the file
    net_handle = Entrez.efetch(
        db="nucleotide", id="NC_006273.2", rettype="fasta", retmode="text"
    )
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
if not os.path.isfile('HCMV.1.bt2'): # To prevent duplication of index
    os.system('bowtie2-build NC_006273.2.fasta HCMV')
os.chdir('..')

# save only the reads that map
os.system(f'bowtie2 -x {p_index}/HCMV -1 {p_data_raw}/SRR5660030_1.fastq.gz -2 {p_data_raw}/SRR5660030_2.fastq.gz -S {p_data_clean}/HCMVmap_30.sam --al-conc-gz {p_data_clean}/SRR5660030_mapped_%.fq.gz')
os.system(f'bowtie2 -x {p_index}/HCMV -1 {p_data_raw}/SRR5660033_1.fastq.gz -2 {p_data_raw}/SRR5660033_2.fastq.gz -S {p_data_clean}/HCMVmap_33.sam --al-conc-gz {p_data_clean}/SRR5660033_mapped_%.fq.gz')
os.system(f'bowtie2 -x {p_index}/HCMV -1 {p_data_raw}/SRR5660044_1.fastq.gz -2 {p_data_raw}/SRR5660044_2.fastq.gz -S {p_data_clean}/HCMVmap_44.sam --al-conc-gz {p_data_clean}/SRR5660044_mapped_%.fq.gz')
os.system(f'bowtie2 -x {p_index}/HCMV -1 {p_data_raw}/SRR5660045_1.fastq.gz -2 {p_data_raw}/SRR5660045_2.fastq.gz -S {p_data_clean}/HCMVmap_45.sam --al-conc-gz {p_data_clean}/SRR5660045_mapped_%.fq.gz')

# write the number of reads that map before and after bowtie2 to log file
my_log = open('PipelineProject.log','w')
for i in accessions:
    # loop through accessions and count read pairs from before and after
    # the counting method is basically dividing number of lines by 4, from here: https://www.biostars.org/p/139006/
    # the subprocess method comes from this thread: https://stackoverflow.com/questions/3503879/assign-output-of-os-system-to-a-variable-and-prevent-it-from-being-displayed-on
    # it allows me to use the command line to count number of reads from each set of files
    # need to decode the output and remove newline character
    bcount = subprocess.check_output(f'echo $(zcat {p_data_raw}/{i}_1.fastq.gz|wc -l)/4|bc', shell=True)
    bcount = bcount.decode('utf-8').strip('\n')
    acount = subprocess.check_output(f'echo $(zcat {p_data_clean}/{i}_mapped_1.fq.gz|wc -l)/4|bc', shell=True)
    acount = acount.decode('utf-8').strip('\n')
    # write to output file
    my_log.write(f'{i} had {bcount} read pairs before Bowtie 2 filtering and {acount} read pairs after\n')
    


# PART 3 ----
#take the bowtie2 output reads and assemble all four genomes
os.chdir(p_data_clean)
# store the individual sub-pairs that will be used for the final SPAdes command
p1 = '--pe-1 1 SRR5660030_mapped_1.fq.gz --pe-2 1 SRR5660030_mapped_2.fq.gz'
p2 = '--pe-1 2 SRR5660033_mapped_1.fq.gz --pe-2 2 SRR5660033_mapped_2.fq.gz'
p3 = '--pe-1 2 SRR5660044_mapped_1.fq.gz --pe-2 2 SRR5660044_mapped_2.fq.gz'
p4 = '--pe-1 2 SRR5660045_mapped_1.fq.gz --pe-2 2 SRR5660045_mapped_2.fq.gz'
# use an f-string to call the correct SPAdes command to create the whole assembly
os.system(f'spades.py -k 77,99,127 -t 2 --only-assembler {p1} {p2} {p3} {p4} -o ../{p_out}/')
# write spades command to output file
my_log.write(f'spades.py -k 77,99,127 -t 2 --only-assembler {p1} {p2} {p3} {p4} -o ../{p_out}/\n')
# my_log.close()
os.chdir('..')

# PART 4 ----
# write python code to calculate the number of contigs > 1000
# Idea for this from here: https://www.biostars.org/p/48797/
# Read in the fasta file
contigs = SeqIO.parse(f'{p_out}/contigs.fasta', 'fasta')
c_filt = []
# If the seq length is greater than 1000, store the seq into a new file contig_filt.fasta
start = 0
for c in contigs:
    if start ==0:
        longest = c
        start +=1
    if len(c.seq) >= 1000:
        c_filt.append(c)
    if len(longest.seq) < len(c.seq):
        longest = c
        
SeqIO.write(c_filt, f'{p_out}/contigs_filt.fasta', 'fasta')
SeqIO.write(longest, f'{p_out}/c_long.fasta', 'fasta')
# Write number of fasta to the log
# my_log = open('PipelineProject.log','w')
scount = subprocess.check_output(f'grep -c "^>" {p_out}/contigs_filt.fasta', shell=True)
scount = scount.decode('utf-8').strip('\n')
# write the number of base pairs to the log
# unix adapted from: https://www.biostars.org/p/78043/#78051
ccount = subprocess.check_output(f"cat {p_out}/contigs_filt.fasta | paste - - - - | cut -f 2 | tr -d '\n' | wc -c", shell=True)
ccount = ccount.decode('utf-8').strip('\n')
# write both to the log file
my_log.write(f'there are {scount} contigs > 1000 bp in the assembly\n')
my_log.write(f'there are {ccount} contigs in the assembly\n')


# PART 5 ----
# Download betaherpesvirinae seqs
# Find entries matching the query
entrezQuery = "txid10357[ORGN:exp]"
os.system(f'esearch -db nucleotide -query "{entrezQuery}" | efetch -format fasta > {p_blast}/blast.fasta')

# Make blast database
os.system(f'makeblastdb -in {p_blast}/blast.fasta -out {p_blast}/BPvirus -title BPvirus -dbtype nucl')

# Make blast query
input_file = f'{p_out}/c_long.fasta'
output_file = f'{p_blast}/myresults.csv'
# using the formatting requested
formatting = '10 sacc pident length qstart qend sstart send bitscore evalue stitle'
os.system(f'blastn -query {input_file} -db {p_blast}/BPvirus -out {output_file} -outfmt "{formatting}"')

# output to the log file
my_log.write(formatting[3:].replace(' ', '\t') + '\n')
my_log.close()
# convert the csv to tsv and append to PipelineProjet.log
os.system(f'head {output_file} |sed "s/,/\t/g" >> PipelineProject.log')





