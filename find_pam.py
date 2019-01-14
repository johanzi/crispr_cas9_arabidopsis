#!/usr/env python


import os
import sys
import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein

'''
Description: Find PAM motifs in a DNA sequence and generate
a genebank file containing the DNA sequence and the PAM 
sequences and the protospacer (20 bp sequence upstream 
of the PAM sequence).

Author: Johan Zicola
Date: Mon Jan 14 20:24:53 CET 2019

'''


# Create parser
parser = argparse.ArgumentParser()

parser.add_argument("-f","--fasta", help="provide one DNA sequence in fasta format", type=str)

parser.add_argument("-o","--output", help="output name, per default name of the fasta file", type=str)

args = parser.parse_args()

# Get name of output
if args.output:
	name_output = os.path.splitext(args.output)[0] + ".gb"
else:
	name_output = os.path.splitext(args.fasta)[0] + ".gb"

# Open argument fasta file
input_handle = open(args.fasta, "rU")

output_handle = open(name_output, "w")

# Check if several fasta sequences in fasta file
sequences = list(SeqIO.parse(input_handle, "fasta"))

if len(sequences) > 1:
	sys.exit("The input fasta file contains "+str(len(sequences))+". I should contain only one sequence")

sequences[0].seq.alphabet = generic_dna

count = SeqIO.write(sequences, output_handle, "genbank")

output_handle.close()
input_handle.close()

print("Coverted %i records" % count)


