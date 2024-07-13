
#Importing Libraries
import re
import os
from Bio import SeqIO
import gzip
import numpy as np
import pandas as pd

#Opening sequence

lines = []

folder_path = '/share/Bacteria/ecoli'
def load_sequences(folder_path):
    sequences = []
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"The folder path '{folder_path}' does not exist.")
    
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.gz'):
            file_path = os.path.join(folder_path, file_name)
            #print(f"Processing file: {file_path}")  # Debugging info
            concatenated_sequence = ''
            with gzip.open(file_path, 'rt') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    concatenated_sequence += str(record.seq)
            sequences.append(concatenated_sequence)
    return sequences

ecoli  = load_sequences(folder_path)
streptococcus = load_sequences('/share/Bacteria/streptococcus')

lines = ecoli + streptococcus
print(lines[0][:50])
print(len(lines))
print(len(list(set(lines))))

#Making reverse compliment
def rev_comp(sequence):
        switch = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
        rev_lines = []

        for i in sequence:
                if i in switch.keys():
                        rev_lines.append(switch[i])
                else:
                        rev_lines.append(i)
        rev_lines = ''.join(rev_lines)
        return rev_lines[::-1]

#Searchers
searcher = "AT.GTA.[A/C]A"
complete_search = "GAGGAA..TC....C.{1,600}?AG....AT.{10,60}?ACA.AA....G..TA"

#Create a list with all of the same sequence
def holder_gen(pssm, sequence1):
        holder = []
        for i in re.finditer(pssm, sequence1):
                holder.append(i.group(0))
        if len(holder) == 0:
                holder.append('-')
        return holder

#Create an array to hold the probabilites
def array_gen(length):
        type = np.array([0, 'A', 'T', 'C', 'G', '-'])
        for i in range(0, length):
                new = np.array([i+1, 0, 0, 0, 0, 0])
                type = np.column_stack((type, new))
        return type

#Add to the array based on each sequence in a list
def pos_matrix(sequences, array):
        for seq in sequences:
                for nuc in range(1, len(seq)+1):
                        index_nuc = dna_to_num(seq[nuc-1])
                        array[index_nuc][nuc] = str(int(array[index_nuc][nuc]) + 1)
        total = len(sequences)
        for i in range(1, array.shape[0]):
                row = array[i][1:array.shape[1]]
                float_array = row.astype(np.float)

                for j in range(0, float_array.shape[0]):
                    array[i][j+1] = format(float_array[j] / total, '.3f')
        return array

#Returns what nucleotide a char is as a number for array placement
#A:1, T:2, C:3, G:4, -:5

def dna_to_num(nuc):
        dna_num = {'A':1, 'T':2, 'C':3, 'G':4, '-':5}
        return dna_num.get(nuc)

#Takes in a  string and calculates how likely it is
def probibility(string, array):
        prob = 1
        for i in range(1, len(string)+1):
                nuc_index = dna_to_num(string[i-1])
                prob *= float(array[nuc_index][i])
        return prob

#Splices off x nucleotides from each string in a list

def cut_off(str_list, num_nuc):
        new_list = []
        for i in range(0, len(str_list)):
                new_list.append(str_list[i][:num_nuc])
                str_list[i] = str_list[i][num_nuc:]
        return new_list, str_list
#Get a list of all reversed and non-revered bacterias
sequences = []
for i in lines:
        sequences.append(i)
        sequences.append(rev_comp(i))
print(len(sequences))


#complete_search = "GAGGAA..TC....C.{1,600}?AG....AT.{10,60}?ACA.AA....G..TA"
#Get Sequences For Forwards
forward_holder = []
for i in range(0, len(sequences)):
        forward_holder.append(holder_gen(complete_search, sequences[i])[0])
print("Initial holder", len(forward_holder))
print("No searcher",  forward_holder.count('-'))

#Remove the strings that don't have a searcher

forward_holder = [i for i in forward_holder if i != '-'] 
 


#Split the sequences so that the {1-600} and {10,60} nucleotide sequences are removed 
#First 16 nuc
fwd_first, middle_holder  = cut_off(forward_holder, 16)


#print("First 16 Nuc", fwd_first)
#Finding the secound part of strand and removing {1, 600} nucs
semi_search = "AG....AT.{10,60}?ACA.AA....G..TA"
semi_holder = []
for i in middle_holder:
        semi_holder.append(holder_gen(semi_search, i)[0])

print("Semi Holder", len(semi_holder))

#Getting next 9 nucs
fwd_sec, middle_holder = cut_off(semi_holder, 9)
#print("Forward Secound",fwd_sec)

#Getting last 15 nucs
fwd_last = []
for i in range(0, len(middle_holder)):
        fwd_last.append(middle_holder[i][-15:])
#print("Forward_last", fwd_last)


#print("GAGGAA..TC....C.{1,600}?AG....AT.{10,60}?ACA.AA....G..TA")
#print(fwd_first + fwd_sec + fwd_last)


#Recombine the sequences with - to indicate nuc gap

whole_sequences = []
for i in range(0, len(fwd_first)):
        whole_sequences.append(str(fwd_first[i] + "-" + fwd_sec[i] + "-" + fwd_last[i]))
#print(whole_sequences)

#Get the matrix array and calculate
array = array_gen(len(whole_sequences[0]))
print(array)
forward_pssm = pos_matrix(whole_sequences, array)


forward_df = pd.DataFrame(forward_pssm)
with pd.option_context('display.max_columns', None,):
    print(forward_df)
