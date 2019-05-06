# Python-Assignment

This repository has two python files.

python_assignment.py : This file contains 7 functions and a "main" function so that a sequence can be analyzed. The script should be run on command line, with the first argument being a filename containing a nucleotide sequence. Example:
$python python_assignment.py test_sequence.txt
The end product of running this script should be:
 (1) the first five lines of a dataframe that contains all possible k values, possible kmers, and observed kmers for the given sequence, 
 (2) the linguistic complexity,
 (3) a pdf file saved into the working directory that has a plot of k vs. proportion of kmers that were observed

test_python_assignment.py : This file contains 8 functions which test the functions in the python_assignment.py file. This file should be used by going to the directory that has this file on the command line and running "pytest".

