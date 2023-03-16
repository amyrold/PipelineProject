#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:27:05 2023

@author: aaronmyrold
"""

import os
from Bio import Entrez
from Bio import SeqIO

# determine the path of the directory this file is located in
# idea taken from here: https://www.pythonanywhere.com/forums/topic/13464/
#my_env = os.path.join(os.path.dirname(__file__)) 
my_env = '/Users/aaronmyrold/Desktop/AM_env/PipelineProject_aaron_myrold'
# set the current working directory to that folder so that remaining paths can function properly
os.chdir(my_env)

# Make directories and paths ----
folder_names = ('1_data_raw', '2_data_clean', '3_output', '4_index', '5_data_test')
p_data_raw = folder_names[0]
p_data_clean = folder_names[1]
p_out = folder_names[2]
p_index = folder_names[3]
p_test = folder_names[4]

# this is primarily due to github not saving empty directories
# need to be able to create folder structure from scratch if needed  
for i in folder_names:
    if not os.path.exists(i):
        os.makedirs(i)


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
    # if switching from test to raw, move test data back to test folder
    os.system(f'mv {p_data_raw}/*.gz {p_test}')
    
    # set working directory to data_raw
    os.chdir('data_raw')
    
    # download all fastq files from ebi database
    # ebi allows you to download the proccess fastq files directly
    os.system('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/000/SRR5660030/SRR5660030_1.fastq.gz')
    os.system('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/000/SRR5660030/SRR5660030_2.fastq.gz')
    
    os.system('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/003/SRR5660033/SRR5660033_1.fastq.gz')
    os.system('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/003/SRR5660033/SRR5660033_2.fastq.gz')
    
    os.system('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/004/SRR5660044/SRR5660044_1.fastq.gz')
    os.system('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/004/SRR5660044/SRR5660044_2.fastq.gz')
    
    os.system('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/005/SRR5660045/SRR5660045_1.fastq.gz')
    os.system('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/005/SRR5660045/SRR5660045_2.fastq.gz')
    
# If we want test data, move it to correct location
if raw_or_test == 'test':
    # move sample data to data_raw folder so that downstream code runs smoothly
    os.system(f'mv {p_test}/*.gz {p_data_raw}')



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
record = SeqIO.read(filename, "fasta") # open file as fasta - i dont think this is needed
os.system('bowtie2-build NC_006273.2.fasta HCMV')
os.chdir('..')

# save only the reads that map
os.system(f'bowtie2 -x {p_index}/HCMV -1 {p_data_raw}/SRR5660030_1.fastq.gz -2 {p_data_raw}/SRR5660030_2.fastq.gz -S {p_data_clean}/HCMVmap_30.sam --al-conc-gz {p_data_clean}/SRR5660030_mapped_%.fq.gz')
os.system(f'bowtie2 -x {p_index}/HCMV -1 {p_data_raw}/SRR5660033_1.fastq.gz -2 {p_data_raw}/SRR5660033_2.fastq.gz -S {p_data_clean}/HCMVmap_33.sam --al-conc-gz {p_data_clean}/SRR5660033_mapped_%.fq.gz')
os.system(f'bowtie2 -x {p_index}/HCMV -1 {p_data_raw}/SRR5660044_1.fastq.gz -2 {p_data_raw}/SRR5660044_2.fastq.gz -S {p_data_clean}/HCMVmap_44.sam --al-conc-gz {p_data_clean}/SRR5660044_mapped_%.fq.gz')
os.system(f'bowtie2 -x {p_index}/HCMV -1 {p_data_raw}/SRR5660045_1.fastq.gz -2 {p_data_raw}/SRR5660045_2.fastq.gz -S {p_data_clean}/HCMVmap_45.sam --al-conc-gz {p_data_clean}/SRR5660045_mapped_%.fq.gz')

# write the number of reads that map before and after bowtie2 to log file
# for file in raw, count reads
# forloop and count from: https://www.biostars.org/p/139006/
os.chdir(p_data_raw)
os.system('for i in `ls *.fastq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done')
# for file in clean, count reads
os.chdir(p_data_raw)
os.system('for i in `ls *.fastq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done')


# PART 3 ----
#take the bowtie2 output reads and assemble all four genomes
# write spades command to output file


# PART 4 ----
# write python code to calculate the number of contigs > 1000
#


# PART 5 ----







