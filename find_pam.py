#!/usr/env python

'''
Description: Find PAM motifs in a DNA sequence and generate
a genebank file containing the DNA sequence and the PAM 
sequences and the protospacer (20 bp sequence upstream 
of the PAM sequence).

Author: Johan Zicola
Date: Mon Jan 14 20:24:53 CET 2019

'''

import os
import sys
import argparse
from Bio import SeqIO, SeqFeature, SeqUtils
from Bio.Alphabet import generic_dna, generic_protein


###################### PARSER #################################

parser = argparse.ArgumentParser()

parser.add_argument("-f","--fasta", help="provide a fasta file containing one DNA sequence only", type=str)

parser.add_argument("-o","--output", help="output name, per default name of the fasta file", type=str)

parser.add_argument("-l","--length-protospacer", help="Length of the protospacer, default = 20 bp", type=int, default=19)

parser.add_argument("-p","--pam", help="Define PAM sequence, default = NGG", type=str, default="NGG")

args = parser.parse_args()

############################################################


#################### FUNCTIONS ############################


# Function to generate a sequence
def create_feature(sequence, name, start, end):
    if str(name) and int(start) and int(end):
        my_feature_location = SeqFeature.FeatureLocation(start, end)
        my_feature = SeqFeature.SeqFeature(my_feature_location,type=name)
        sequence.features.append(my_feature)

# example: create_feature(sequence, "glup", 20,40)


# Function to find PAM
def search_motif(sequence):

    motif = str(argparse.pam)
    
    len_motif = int(len(motif))

    len_protospacer = int(argparse.length_protospacer)
    
    full_len = len_motif + len_protospacer

    len_dna = int(len(sequence.seq))

    # Output of nt_search is a list containing the motif and the start position (0-based)
    # of every hit in the DNA sequence

    # Search on fw strand
    matches_fw = SeqUtils.nt_search(str(sequence.seq), motif)
    
    start_positions_fw = matches_fw[1::]

    end_positions_fw = start - full_len for start in start_positions_fw

    # Loop over positions and create end position and generate a feature
    # based on length of protospacer


    # The coordinates are different and need to be corrected to match to fw strand
    matches_rv = SeqUtils.nt_search(str(sequence.seq.reverse_complement()), motif)
    
    start_positions_rv = matches_fw[1::]

    end_positions_rv = start - full_len for start in start_positions_rv

    # Need to convert the coordinates in forward strand
    start_positions_rv = len_dna - start for start in start_positions_rv
    end_positions_rv = len_dna - end for end in end_positions_rv

    # Either I create all features here or I just provide a tuple of 4 lists 
    return 




############################################################


# Get name of output
if args.output:
	name_output = os.path.splitext(args.output)[0] + ".gb"
else:
	name_output = os.path.splitext(args.fasta)[0] + ".gb"


# Open argument fasta file
input_handle = open(args.fasta, "rU")

# Open output file to write in
output_handle = open(name_output, "w")


# Check if several fasta sequences in fasta file
sequences = list(SeqIO.parse(input_handle, "fasta"))

if len(sequences) > 1:
	sys.exit("The input fasta file contains "+str(len(sequences))+". I should contain only one sequence")
else:
    sequence = sequences[0]


# Get DNA sequence from the fasta file
sequence.seq.alphabet = generic_dna


# Find PAM in DNA sequence
# https://biopython.org/wiki/Seq

search_motif(sequence)



# Add features before converting to genbank

# Convert fasta file into genebank file. Not
# Note that conversion to genbank format leads to the addition of 1 for 
# each start of features. Therefore, the start position should be given 
# 0-based system to anticipate the 1-based conversion
SeqIO.write(sequence, output_handle, "genbank")



output_handle.close()
input_handle.close()

# Get positions of the PAM and corresponding protospacers
# The genbank format is 1-based, closed [start, end], range = end - start + 1







