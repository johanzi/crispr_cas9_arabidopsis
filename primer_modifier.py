
# Libaries
import sys

# Usage
"""
CRISPR-Cas9 cloning based on paper Wang et al 2015 (Genome Biology)

Generate primers with indicated target sequence for BsaI and overhang,
target sites, and other primers with overlap with pCBC-DT1T2 vector
The output is the generation of 4 primer sequences to be used for 
PCR reaction with pCBC-DT1T2 for the 2 target sites which can be incorporated

"""


def reverse_complement(my_dna):
    my_dna = str(my_dna)
    my_dna = my_dna.upper()
    new_str1 = str.replace(my_dna, 'A','t')
    new_str2 = str.replace(new_str1, 'T','a')
    new_str3 = str.replace(new_str2, 'C','g')
    new_str4 = str.replace(new_str3, 'G','c')
    new_str5 = new_str4.lower()
    new_str6 = new_str5[::-1]

    return new_str6


# Define sequence to add for sgRNA1 and 2
# User needs to fill in according to requested overhang
# for Golden Gate cloning and used pCBC vector


length_protospacer = 20

# sgRNA1 (first target site)
golden_gate_1_5p = "ATATATGGTCTCGATTG"
golden_gate_1_3p = "G"

vector_overhang_1 = "GTTTTAGAGCTAGAAATAGCA"

# sgRNA2 (second target site)
# Note that in Wang et al 2015, the primer DT2-BsR contain an error
# of 1 base pair ATTATTGGTCTC|G|AAAC where G should be T. It actually
# should not influence the construct since this part of the sequence is 
# cut by BsaI
golden_gate_2_5p = "ATTATTGGTCTCTAAAC"
golden_gate_2_3p = "CAA"

vector_overhang_2_5p = "AAC"
vector_overhang_2_3p = "CAATCTCTTAGTCGACTCTAC"

# Get input file containing primer sequence
name_file = sys.argv[1] 

file = open(name_file, "r")


# Add prefix to primer name

if len(sys.argv) == 3:
    prefix = str(sys.argv[2])+"_"
else:
    prefix = ""


# Initialize sequence number
num_sequence = 1

for line in file:
    if line.strip():
        line = line.strip()
        
        if len(line) != length_protospacer:
            sys.exit("The target site "+str(line) +"  is not equal to "+str(length_protospacer))
        else:
            line = line.lower()
            
            DT1_BsF = golden_gate_1_5p + str(line) + golden_gate_1_3p
            DT1_F0 = str(line) + vector_overhang_1
            
            line_rc = reverse_complement(line)
			
            DT2_BsR = golden_gate_2_5p + str(line_rc) + golden_gate_2_3p
            DT2_R0 = vector_overhang_2_5p + str(line_rc) + vector_overhang_2_3p

            print("Target sequence "+str(num_sequence)+": "+line)
            print(prefix+"DT1_BsF\t"+DT1_BsF)
            print(prefix+"DT1_F0\t"+DT1_F0)
            print(prefix+"DT2_BsR\t"+DT2_BsR)
            print(prefix+"DT2_R0\t"+DT2_R0)
            print("\n")

            num_sequence += 1

