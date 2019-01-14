#!/usr/env python


import os
import sys
import argparse
from Bio import SeqIO, SeqFeature
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
else:
    sequence = sequences[0]

# Get DNA sequence from the fasta file
sequence.seq.alphabet = generic_dna

# Add features before converting to genbank

# Function to generate a sequence
def create_feature(sequence, name, start, end):
    if str(name) and int(start) and int(end):
        my_feature_location = SeqFeature.FeatureLocation(start, end)
        my_feature = SeqFeature.SeqFeature(my_feature_location,type=name)
        sequence.features.append(my_feature)

create_feature(sequence, "flop", 10,20)
create_feature(sequence, "jesus", 100,200)
create_feature(sequence, "glup", 20,40)


print(sequence.features)

# Convert fasta file into genebank file. Not
# Note that conversion to genbank format leads to the addition of 1 for 
# each start of features. Therefore, the start position should be given 
# 0-based system to anticipate the 1-based conversion
SeqIO.write(sequences, output_handle, "genbank")


# Find PAM in DNA sequence
# https://biopython.org/wiki/Seq
print(sequences[0].seq.find("ACAAGT"))

output_handle.close()
input_handle.close()

# Get positions of the PAM and corresponding protospacers
# The genbank format is 1-based, closed [start, end], range = end - start + 1







