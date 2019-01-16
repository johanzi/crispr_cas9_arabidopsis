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

# Note I need to provide both const and default argument to have these values even if flag is not given
parser.add_argument("-l","--length-protospacer", help="Length of the protospacer, default = 20 bp", nargs='?', type=int, const=19, default=20)

parser.add_argument("-p","--pam", help="Define PAM sequence, default = NGG", nargs='?', type=str, const="NGG", default="NGG" )

args = parser.parse_args()

############################################################



#################### FUNCTIONS ############################


# Function to generate a sequence
# Use option strand to indicate on which strand is located the feature (fw per default)
# Use quote +1 (forward) or -1 (reverse)
# Check http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html

def create_feature(sequence, name, start, end, strand=+1):

    if str(name) and int(start) and int(end):
        my_feature_location = SeqFeature.FeatureLocation(start, end, strand=strand)
        my_feature = SeqFeature.SeqFeature(my_feature_location,type=name)
        sequence.features.append(my_feature)

# example: create_feature(sequence, "glup", 20,40)



# Function to find PAM
# The coordinates given here are 0-based (the conversion to Genbank format
# will add 1bp to the start coordinate and make it correctly matching the sequence)
def search_motif(sequence):

    motif = str(args.pam)
    
    len_motif = int(len(motif))

    len_protospacer = int(args.length_protospacer)
    
    full_len = len_motif + len_protospacer

    len_dna = int(len(sequence.seq))

    # Output of nt_search is a list containing the motif and the start position (0-based)
    # of every hit in the DNA sequence

    # Search on fw strand
    matches_fw = SeqUtils.nt_search(str(sequence.seq), motif)
    
    # Initialyze final list
    coordinates_fw = []
     
    if len(matches_fw) > 1:
        end_positions_fw = matches_fw[1::]
        start_positions_fw = [ end - len_protospacer for end in end_positions_fw ]
        
        # Check if protospacer fits in the sequence before adding the start
        # and end coordinate to the list
        for start, end in zip(start_positions_fw, end_positions_fw):
            if start > 0:
                coordinates_fw.append([start, end])    

    # The coordinates are different and need to be corrected to match to fw strand
    reverse_seq = str(sequence.seq.reverse_complement())
    
    matches_rv = SeqUtils.nt_search(reverse_seq, motif)
    
    # Initialyze final list
    coordinates_rv = []

    if len(matches_rv) > 1:
        end_positions_rv = matches_rv[1::]
        start_positions_rv = [ end - len_protospacer for end in end_positions_rv ]
        # Need to convert the coordinates in forward strand
        end_positions = [ len_dna - start for start in start_positions_rv ]
        start_positions = [ len_dna - end for end in end_positions_rv ]
        
        # Check if protospacer fits in the sequence before adding the start
        # and end coordinate to the list
        for start, end in zip(start_positions, end_positions):
            if start > 0 and end < len_dna:
                coordinates_rv.append([start, end])    
        
    
    # Return a tuple of lists for fw and rv matches
    return coordinates_fw, coordinates_rv


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
coordinates_fw, coordinates_rv = search_motif(sequence)

print(coordinates_fw)
print(coordinates_rv)

# Create feature for each protospacers
if len(coordinates_fw) > 0:
    for index, item in enumerate(coordinates_fw):
        protospacer_name = "site_fw_" + str(index + 1)
        create_feature(sequence, protospacer_name, item[0], item[1], +1)

if len(coordinates_rv) > 0:
    for index, item in enumerate(coordinates_rv):
        protospacer_name = "site_rv_" + str(index + 1)
        create_feature(sequence, protospacer_name, item[0], item[1], -1)



# Convert fasta file into genebank file. Not
# Note that conversion to genbank format leads to the addition of 1 for 
# each start of features. Therefore, the start position should be given 
# 0-based system to anticipate the 1-based conversion
SeqIO.write(sequence, output_handle, "genbank")



output_handle.close()
input_handle.close()

# Get positions of the PAM and corresponding protospacers
# The genbank format is 1-based, closed [start, end], range = end - start + 1







