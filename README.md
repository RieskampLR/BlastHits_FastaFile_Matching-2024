README

Program description:
Automated insertion of hitDescriptions from a blastx file into the headers of a fasta file.
The result is stored in an output file (text format).
Sequences without a blast hit are excluded.
Lea Rachel Rieskamp
14.10.2024

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
