 #versione di test

import sys
import os
import Bio
import pandas as pd
from Bio import SeqIO

#checks for arguments
if len(sys.argv) <= 1:
    print("Please provide FASTA file")

#help message
if sys.argv[1] == '-h' or sys.argv[1] == 'help' or sys.argv[1] == '-help':
    print("Usage: fasta2csv inputfile.fst outputfile.csv")

#INPUT
input = sys.argv[1]

#verifies if the input file exists
if not os.path.exists(input):
    print('Error: file %s does not exist!', input)

#verifies if the output was passed in the arguments
output = 'output.csv'
if len(sys.argv) > 2:
    output = sys.argv[2]

#dataframe definition   
df = pd.DataFrame()

#conversion from FASTA to dataframe via SeqIO
for seq_record in SeqIO.parse(input, "fasta"):
    df = df.append({0: seq_record.id,1:str(seq_record.seq)}, ignore_index = True)
    
#splitting the sequences strings into a second dataframe
M = df[1].str.split("", expand = True)

#dictionary definition
dictionary = {
        'A':27,
        'B':1,
        'C':2,
        'D':3,
        'E':4,
        'F':5,
        'G':6,
        'H':7,
        'I':8,
        'J':9,
        'K':10,
        'L':11,
        'M':12,
        'N':13,
        'O':14,
        'P':15,
        'Q':16,
        'R':17,
        'S':18,
        'T':19,
        'U':20,
        'V':21,
        'W':22,
        'X':23,
        'Y':24,
        'Z':25,
        '*':26,
        '-':0,
}

#applying the dictionary to the dataframe
M = M.replace(dictionary)

#writing the csv
M.to_csv(output)