
# Libaries
import sys

# Usage
"""
CRISPR-Cas9 cloning based on paper Wang et al 2015 (Genome Biology)

Generate primers with indicated target sequence for BsaI and overhang,
target sites, and other primers with overlap with pCBC-DT1T2 vector
The output is the generation of 2 primer sequences to be used for 
PCR reaction with pCBC-DT1T2

"""

# Define sequence to add for sgRNA1 and 2
# User needs to fill in according to requested overhang
# for Golden Gate cloning and used pCBC vector

# sgRNA1
golden_gate_1_5p = "ATATATGGTCTCGATTG"
golden_gate_1_3p = "G"

vector_overhang_1 = "GTTTTAGAGCTAGAAATAGCA"

# sgRNA2
golden_gate_2_5p = "ATTATTGGTCTCTAAAC"
golden_gate_2_3p = "CAA"

vector_overhang_2_5p = "AAC"
vector_overhang_2_3p = "CAATCTCTTAGTCGACTCTAC"

# Get input to know if target site is in sgRNA1 or sgRNA2
target_site = sys.argv[1] 


# Get input file containing primer sequence
name_file = sys.argv[2] 

file = open(name_file, "r")

for line in file:
    line = line.strip()
    
    if len(line) != 19:
        sys.exit("The target site "+str(line) +"  is not equal to 19")
    else:
        if target_site == "1":
            primer_golden_gate = golden_gate_1_5p + str(line) + golden_gate_1_3p
            primer_vector = str(line) + vector_overhang_1
        elif target_site == "2":
            primer_golden_gate = golden_gate_2_5p + str(line) + golden_gate_2_3p
            primer_vector = vector_overhang_2_5p + str(line) + vector_overhang_2_3p
        else:
            sys.exit("Choose either '1' or '2' for target site 1 and 2, respectively")

        print("Target site: "+line)
        print("sgRNA site: "+str(target_site))
        print("golden gate primer: "+ primer_golden_gate)
        print("vector primer: "+primer_vector)
