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
my_env = os.path.join(os.path.dirname(__file__)) 
#my_env = '/Users/aaronmyrold/Desktop/PipelineProject/PipelineProject_aaron_myrold'
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
        
# Make accession list ----
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
    os.system('fasterq-dump --threads 2 --progress SRR5660030')
    os.system('fasterq-dump --threads 2 --progress SRR5660033')
    os.system('fasterq-dump --threads 2 --progress SRR5660044')
    os.system('fasterq-dump --threads 2 --progress SRR5660045')
    
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
# for file in raw, count reads

my_log = open('PipelineProject.log','w')
# bcount = {}
# acount = {}
for i in accessions:
    bcount = os.system(f'echo $(zcat {p_data_raw}/{i}_1.fastq.gz|wc -l)/4|bc') + os.system(f'echo $(zcat {p_data_raw}/{i}_2.fastq.gz|wc -l)/4|bc')
    acount = os.system(f'echo $(zcat {p_data_clean}/{i}_mapped_1.fq.gz|wc -l)/4|bc') + os.system(f'echo $(zcat {p_data_clean}/{i}_mapped_2.fq.gz|wc -l)/4|bc')
    my_log.write(f'{i} had {bcount} read pairs before Bowtie 2 filtering and {acount} read pairs after')
    
    


# # forloop and count from: https://www.biostars.org/p/139006/
# os.chdir(p_data_raw)
# os.system('for i in `ls *.fastq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done')
# # for file in clean, count reads
# os.chdir(f'../{p_data_clean}')
# os.system('for i in `ls *.fq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done')


# # PART 3 ----
# #take the bowtie2 output reads and assemble all four genomes
# # write spades command to output file
# os.chdir(p_data_clean)
# p1 = '--pe-1 1 SRR5660030_mapped_1.fq.gz --pe-2 1 SRR5660030_mapped_2.fq.gz'
# p2 = '--pe-1 2 SRR5660033_mapped_1.fq.gz --pe-2 2 SRR5660033_mapped_2.fq.gz'
# p3 = '--pe-1 2 SRR5660044_mapped_1.fq.gz --pe-2 2 SRR5660044_mapped_2.fq.gz'
# p4 = '--pe-1 2 SRR5660045_mapped_1.fq.gz --pe-2 2 SRR5660045_mapped_2.fq.gz'

# os.system(f'spades.py -k 77,99,127 -t 2 --only-assembler {p1} {p2} {p3} {p4} -o ../{p_out}/')




# PART 4 ----
# write python code to calculate the number of contigs > 1000
#


# PART 5 ----








