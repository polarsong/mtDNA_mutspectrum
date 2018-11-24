#!/usr/bin/env python

from Bio import SeqIO
import sys
from collections import Counter

"""
	This script parses the provided .gb source file and creates a file of features with header string "species	taxonomy	A	T	G	C	X	FULL_GENOME" named full_table.txt unless otherwise specified

	Usage: generate_full_table.py input_genbank_file [name_of_output_file]
"""

# If no input file specified, throw an error
elif len(sys.argv) < 2:
	sys.exit("No input file! Please provide genbank file name as input. Check the docstring for more info.")

user_input = sys.argv[1]
input_file = open(user_input, "rt")


# Check if different output file name specified by user
if len(sys.argv) > 2:
	user_output = sys.argv[2]
	output_file = open(user_output, "wt")
else:
	output_file = open("full_table.txt", "wt")


for rec in SeqIO.parse(input_file, 'genbank'):
	species = "_".join(rec.annotations['organism'].split())
	taxonomy = "_".join(rec.annotations["taxonomy"])
	count = Counter(rec.seq)
	(A, C, G, T) = (count["A"], count["C"], count["G"], count["T"])
	X = str(len(rec.seq) - (A + C + G + T))
	output_file.write("\t".join([species, taxonomy, str(A), str(C), str(G), str(T), X, str(rec.seq)]) + "\n")

input_file.close()
output_file.close()