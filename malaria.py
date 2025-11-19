#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
malaria.py

Description:
This program inserts hitDescriptions from a blastx file
into the header of a fasta file and stores the result in an output file (text format).
Sequences without a blast hit are excluded.

User-defined functions: None
Non-standard modules: None

Procedure:
    1. Input checks
        1.1 Number of command line arguments
        1.2 Existence of files
        1.3 Size of files
        1.4 Files types by name format
        1.5 File types by content characteristics
    2. Iteration over the blastx file to create a dictionary of the hitDescriptions
       using the sequence IDs as keys. 'null'-hits are excluded.
    3. Iteration over the fasta file to write the entries into the output file.
       Only entries which´s ID is found in the dictionary are included.
       The associated hitDescription from the dictionary is added to the entry´s header.

Input: blastx file , fasta file , name of output file
Output: output file

Usage: malaria.py yourfastafile.fna yourblastxfile.blastx.tab youroutputfile.txt

Version: 1.00
Date: 2024-10-14
Author: Lea Rachel Rieskamp

"""

# Imports:

import sys
import os
from pathlib import Path



# Input file check - number of arguments:

if len(sys.argv) != 4:
    print("Incorrect number of command line arguments.\nTry: python malaria.py yourfastafile yourblastxfile youroutputfile\nProgram terminated.")
    exit()



# Storing Paths of arguments in variables:
    
closedoutfile = Path(sys.argv[3])
closedblastfile = Path(sys.argv[2])
closedfastafile = Path(sys.argv[1])



# Input file check - existence:

if closedoutfile.is_file():
    bresponseoutfile = input(f"{closedoutfile} already exists. Continue and overwrite existing file?")
    if not bresponseoutfile:
        bresponseoutfile="empty"
    if bresponseoutfile.lower()[0] == "n":
        exit()

if not closedblastfile.is_file():
    print(f"{closedblastfile} not found.\nProgram terminated.")
    exit()

if not closedfastafile.is_file():
    print(f"{closedfastafile} not found.\nProgram terminated.")
    exit()



# Input file check - size:
# 1GB warning

blastsize = os.path.getsize(closedblastfile)
fnasize = os.path.getsize(closedfastafile)

bigfilesize = 1073741824

if blastsize > bigfilesize:
    bresponsebsize = input(f"{closedblastfile} > 1GB.\nDo you want to continue with this file?")
    if bresponsebsize.lower()[0] == "n":
       exit()

if fnasize > bigfilesize:
    bresponsefsize = input(f"{closedfastafile} > 1GB.\nDo you want to continue with this file?")
    if bresponsefsize.lower()[0] == "n":
       exit()



# Input file check - file types by name format:

if not str(closedblastfile).endswith(".blastx.tab"):
    print(f"{closedblastfile} is not a blastx.tab file. File has to end on .blastx.tab.\nCheck if file is compressed.\nCheck if order of commandline arguments is correct (python malaria.py yourfastafile yourblastxfile youroutputfile)\nProgram terminated.")
    exit()

if not str(closedfastafile).endswith((".fna",".fasta",".fa")):
    if str(closedfastafile).endswith(".fastq"):
        print("You provided a FASTQ file. Please provide a FASTA file instead.\nProgram terminated.")
        exit()
    print(f"{closedfastafile} is not a FASTA file. File has to end on .fna, .fasta, or .fa.\nCheck if file is compressed.\nCheck if order of commandline arguments is correct (python malaria.py yourfastafile yourblastxfile youroutputfile)\nProgram terminated.")
    exit()



# Open output and blast file variables:
    
outfile = open(closedoutfile,"w")
blastfile = open(closedblastfile,"r")



# Input file check - file content:
# blast file start: "#"
# fasta file start: ">"

blastfirstline = blastfile.readline()
if blastfirstline[0] != "#":
    print(f"{closedblastfile} is not a blastx file or has incorrect header.\nProgram terminated.")
    exit()
blastfile.seek(0)
# .seek(0) sets read position back to first line

ffastafile = open(closedfastafile,"r")
fastafirstline = ffastafile.readline()
if fastafirstline[0] != ">":
    ffastafile.close()
    print(f"{closedfastafile} is not a FASTA file. File has to start with '>'.\nProgram terminated.")
    exit()
ffastafile.close()



# Dictionary of IDs(key) with hitDescription(item) from blastxfile:

protdes = {}
    
firstline = 1
for bline in blastfile:
    fields = bline.split("\t")
    if firstline:
        try:
            nHitdescolumn = fields.index("hitDescription")
        except ValueError:
            print(f"Header or Column 'hitDescription' missing in {closedblastfile}.\nProgram terminated")
            sys.exit()
        firstline = 0
    else:
        if fields[nHitdescolumn] != "null":
            protdes[fields[0]] = fields[nHitdescolumn]

# Explanation:
# 1. Creates empty dictionary
# 2. Sets firstline true
# 3. Loops through blastfile line by line
# 4. First line will go into try-except-block,
#    tests if hitDescriptions present
#    -if yes, column number stored in nHitdescolumn, and firstline set false
#    -if not, error message given and program stopped
# 5. Next lines skip the if-firstline-block
#    -if hitDescription found in the blast search (hitDescription is not null),
#     it is stored in the dictionary with its associated ID as the key
#    -if not, loop continues with next line



# Writing headers, hitDescriptions, and sequences into output file:

count = 0
totalcount = 0 

with open(closedfastafile,"r") as fastafile:
    current_header = ''
    seq = ''
    for fline in fastafile:
        if ">" in fline:
            totalcount += 1
            if current_header:
                fsplitline = current_header.split("\t")
                seqid = fsplitline[0][1:]
                if seqid in protdes:
                    hline = f"{current_header.strip()}\tprotein={protdes[seqid]}\n{seq}\n"
                    ## print(hline)
                    count += 1
                    outfile.write(hline)
            current_header = fline
            seq = ''

        else:
            seq += fline.strip()
            
    fsline = current_header.split("\t")
    seqid = fsline[0][1:]
    if seqid in protdes:
        hline = f"{current_header.strip()}\tprotein={protdes[seqid]}\n{seq}\n"
        count += 1
        outfile.write(hline)

# Explanation:
# 1. opens file and loops line by line
# 2. On encounter of first header (">"), goes into first if block,
#    current_header is empty, so skips second if block,
#    and directly assigns the first header into the variable current_header
# 3. Next line, a sequence, does not have ">", hence skips if block, and goes into else block,
#    where the sequence is stored in the variable seq
# 4. Next line
#    - if multi-line sequence: The line contains the sequence continuation,
#      no ">", so skips the first if block
#      goes to else block and is added to the first part of the sequence (seq +=)
#    - if one-line sequence: The line contains the next header,
#      has ">" and goes into first if block
#      current_header contains the first header, so not empty, so goes into second if block
#      - if the first header ID is in the protdes disctionary (so not a 'null' hit),
#        the first header, the belonging protein description, and the sequence are written into the output file
#        then the current_header variable is updated with the next header
#        and the seq variable is emptied
#      - if the first header is not in the protdes dictionary
#        the current_header variable is directly updated with the next
#        and the seq variable emptied
# 5. Repeats line by line
# 6. Reaching the last line, the variables are still updated and stored but the loop is not reentered,
#    so when loop is exited, the variables are written separately into the output file,
#    if the header ID is in the protdes disctionary
      


# Closing open files:
    
blastfile.close()
outfile.close()
    





    
    
    

