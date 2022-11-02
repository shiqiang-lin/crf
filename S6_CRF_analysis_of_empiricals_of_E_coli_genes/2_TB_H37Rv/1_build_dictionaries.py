"""
The py file creates the dictionaries for CRF calculation, which includes

(1) codon_to_amino_acid_dictionary
(2) amino_acid_to_codons_dictionary
(3) state_pair_dictionary
(4) state_observation_pair_dictionary
(5) observation_states_dictionary

The dictionaries (3) and (4) actually describe the graph structure of the CRF model.

To run the build_dictionaries.py, the codon_table_raw.txt file, which has been manually built
and stores the codon information, is required.
"""
import numpy as np


#get the codon file
codon_txt_file = open('0_codon_table_raw.txt','r')
codon_txt_file_lines = codon_txt_file.readlines()
codon_txt_file.close()

#print to take a look
print("\n")
print("Here is the codon table:\n")

for i in range(len(codon_txt_file_lines)):
    codon_txt_file_line = codon_txt_file_lines[i].strip()
    print(codon_txt_file_line)

print("\n")

###dictionary (1) the codon to amino acid table for gene translation
codon_to_amino_acid_dictionary = {}

for i in range(len(codon_txt_file_lines)):
    codon_txt_file_line = codon_txt_file_lines[i].strip()
    codon_txt_file_line_split = codon_txt_file_line.split()

    amino_acid = codon_txt_file_line_split[0]
    codon_number = len(codon_txt_file_line_split) - 1

    for i in range(codon_number):
        codon_to_amino_acid_dictionary[codon_txt_file_line_split[i+1]] = amino_acid

#show the codons for gene translation
print("\n")
print("Here is the codon_to_amino_acid_dictionary.")
#print(codon_to_amino_acid_dictionary)
"""
for item in codon_to_amino_acid_dictionary.items():
    print(item)
"""
print(codon_to_amino_acid_dictionary)
print("There are %s codons in total." % len(codon_to_amino_acid_dictionary))
print("\n")


#save dictionary (1) the codon to amino acid table. #The .npy will be added to the filename.
np.save('codon_to_amino_acid.dict.npy',codon_to_amino_acid_dictionary)                         




"""       
For constructing dictionary (3), we need first to list all the amino acid pairs and then find all the corresponding
codon combinations. To do this, we store the txt lines in the dictionary for better accessibility. This dictionary
is also useful for forward and backward algorithms later.
"""
###(2) the amino acid to codons dictionary
amino_acid_to_codons_dictionary = {}

for i in range(len(codon_txt_file_lines)):
    codon_txt_file_line = codon_txt_file_lines[i].strip()
    codon_txt_file_line_split = codon_txt_file_line.split()
    amino_acid_to_codons_dictionary[codon_txt_file_line_split[0]] = codon_txt_file_line_split[1:]
                                                                                
print("Here is the amino_acid_to_codons_dictionary.")
print(amino_acid_to_codons_dictionary)
print("\n")


#save dictionary (2) the amino acid to codons dictionary. #The .npy will be added to the filename.
np.save('amino_acid_to_codons.dict',amino_acid_to_codons_dictionary)          



"""
******************************************__dictionary (3)__begin__****************************************************
Here we use the nested dictionary to show the amino acid pair - state pairs. Also, the print code is included to see
the resultant nested dictionary. The nested dictionary is long and we can copy to TxtFile to take a look.

The observation pairs of the dictionary consist of 
1)observation_0 observation_1, in which the observation_0 is START and the observation_1 is a start amino acid;
2)observation_1 observation_2, in which the observation_2 is a normal amino acid;
3)observation_i observation_i+1 (1<i<n), in which observation_i is a normal amino acid and the observation_i+1 is a
  normal amino acid or a *(stop);
4)observation_n observation_n+1, in which the observation_n is a * and the observation_n+1 is STOP.

Of these parts, the 3) and 4) are invariable, however, 1) and 2) depend on the codons of first amino acids of the
proteins encoded by the genes of the genome.
"""
##This part deals with 3) of the dictionary.
"""
There are 20 lefts. The '*' only stays in the end of protein sequence.
"""
lefts_of_amino_acid_pair = []

for i in range(len(codon_txt_file_lines)-1):
    codon_txt_file_line = codon_txt_file_lines[i].strip()
    codon_txt_file_line_split = codon_txt_file_line.split()
    lefts_of_amino_acid_pair.append(codon_txt_file_line_split[0])

print("Here are the lefts of amino acid pair.")
print(lefts_of_amino_acid_pair)
print("There are %s lefts of amino acid pair in total." %len(lefts_of_amino_acid_pair))
print("\n")

""" 
There are 21 rights, including '*'.
"""
rights_of_amino_acid_pair = []

for i in range(len(codon_txt_file_lines)):
    codon_txt_file_line = codon_txt_file_lines[i].strip()
    codon_txt_file_line_split = codon_txt_file_line.split()
    rights_of_amino_acid_pair.append(codon_txt_file_line_split[0])
    
    
print("Here are the rights of amino acid pair.")
print(rights_of_amino_acid_pair)
print("There are %s rights of amino acid pair in total." %len(rights_of_amino_acid_pair))
print("\n")

amino_acid_pairs_normal = []

for i in lefts_of_amino_acid_pair:
    for j in rights_of_amino_acid_pair:
        amino_acid_pairs_normal.append(i+j)

print("Here are the amino acid pairs of the amino_acid_pair_normal.")
print(amino_acid_pairs_normal)
print("There are %s pairs of amino_acid_pair_normal in total." %len(amino_acid_pairs_normal))
print("\n")


state_pair_dictionary = {}          #outer dictionary
codon_pair_dictionary = {}          #inner dictionary

for left_right in amino_acid_pairs_normal:
    left = left_right[0]
    right = left_right[1]

    left_codons_list = amino_acid_to_codons_dictionary[left]
    right_codons_list = amino_acid_to_codons_dictionary[right]
  

    for left_codon in left_codons_list:
        for right_codon in right_codons_list:
            codon_pair = left_codon + right_codon
            #print(codon_pair)    
            codon_pair_dictionary[codon_pair] = [0, 0, 0]

    state_pair_dictionary[left_right] = codon_pair_dictionary
    codon_pair_dictionary = {}      #need to empty the inner dictionary now


##This part deals with 4) of the dictionary.
state_pair_dictionary['*STOP'] = {'TAASTOP':[0, 0, 0], 'TAGSTOP':[0, 0, 0], 'TGASTOP':[0, 0, 0]}


"""
Now deal with the 1) and 2), which might be varible with different genomes. The start codons of all genes shall be
checked before running the program. This is done in the gene preparation .py file. The start codons can be added
to the dictionary if rare situations occur.
"""
#This part deals with 1) of the dictionary.
state_pair_dictionary['STARTM'] = {'STARTATG':[0, 0, 0]}
state_pair_dictionary['STARTVs'] = {'STARTGTG':[0, 0, 0]}
state_pair_dictionary['STARTLs'] = {'STARTTTG':[0, 0, 0], 'STARTCTG':[0, 0, 0]}
state_pair_dictionary['STARTIs'] = {'STARTATT':[0, 0, 0]}


#This part deals with 2) of the dictionary.
"""
add the state pairs consisting of start codons Vs, Ls and Is. This part may be modified according to different
genomes with genes of various start codons.
"""
first_amino_acids = ['Vs','Ls','Is']
amino_acid_pairs_special = []

for i in first_amino_acids:
    for j in lefts_of_amino_acid_pair:
        amino_acid_pairs_special.append(i+j)

print("Here are amino acid pairs of the amino_acid_pairs_special.")
print(amino_acid_pairs_special)
print("There are %s pairs of amino_acid_pairs_special in total." %len(amino_acid_pairs_special))
print("\n")


first_amino_acid_to_codons_dictionary = {'Vs':['GTG'], 'Ls':['TTG', 'CTG'], 'Is':['ATT']}

for left_right in amino_acid_pairs_special:
    left = left_right[0:2]
    right = left_right[2]

    left_codons_list = first_amino_acid_to_codons_dictionary[left]
    right_codons_list = amino_acid_to_codons_dictionary[right]
  

    for left_codon in left_codons_list:
        for right_codon in right_codons_list:
            codon_pair = left_codon + right_codon
            #print(codon_pair)    
            codon_pair_dictionary[codon_pair] = [0, 0, 0]

    state_pair_dictionary[left_right] = codon_pair_dictionary
    codon_pair_dictionary = {}      #need to empty the inner dictionary now



#print to see the state_pair_dictionary
print("Here is the state_pair_dictionary.")
print(state_pair_dictionary)
print("The length of the state_pair_dictionary is %s." %len(state_pair_dictionary))
print("\n")




#save dictionary (3) state_pair_dictionary. #The .npy will be added to the filename.
np.save('state_pair.dict',state_pair_dictionary) 

"""
**********************************************__dictionary__(3)__end__*************************************************
"""



"""
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$__dictionary (4)__begin__$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
dictionary (4) the state-observation pair (codon amino acid) dictionary
The observations of the dictionary (4) consist of
1) normal amino acids;
2) start amino acids;
3) START and STOP.
"""
##This part deals with 1) of the dictionary.
state_observation_pair_dictionary = {}

for i in range(len(codon_txt_file_lines)):
    codon_txt_file_line = codon_txt_file_lines[i].strip()
    codon_txt_file_line_split = codon_txt_file_line.split()

    states = codon_txt_file_line_split[1:]
    observation = codon_txt_file_line_split[0]

    for state in states:
        state_observation_pair_dictionary[state+observation] = [0,0,0]
        

"""
The first amino acid is considered as different from other amino acids in the middle of the sequence if the number
or types of the states are not the same. Therefore, they are manually added to the state_observation_pair_dictionary.

If not differentiated, when calculating the empirical features, it is not possible to make difference
codoni-AA1(first codon and first amino acid) with codoni-AAj(1<j<n). Then the parameters of the codoni-AA1 will be
affected by the codoni-AAj. To solve this problem, the AA1 is regarded as a new observation if the number or types of
the states are not the same with that in the middle of the sequence. 
"""

##This part deals with 2) of the dictionary.
"""
The 's' means start amino acid. 's' in state_observation_pair_dictionary['ATGMs'] = [0, 0, 0] is not added,
because M has only one state, i.e., ATG, which is the same as M in the middle of the sequence.
"""
state_observation_pair_dictionary['ATGM'] = [0,0,0]
state_observation_pair_dictionary['GTGVs'] = [0,0,0]
state_observation_pair_dictionary['TTGLs'] = [0,0,0]
state_observation_pair_dictionary['CTGLs'] = [0,0,0]
state_observation_pair_dictionary['ATTIs'] = [0,0,0]

##This part deals with 3) of the dictionary.
state_observation_pair_dictionary['STARTSTART'] = [0,0,0]
state_observation_pair_dictionary['STOPSTOP'] = [0,0,0]



#print to see the state_observation_pair_dictionary
print("\n")
print("Here is the state_observation_pair_dictionary.")
print(state_observation_pair_dictionary)
print("The length of the state_observation_pair_dictionary is %s." % len(state_observation_pair_dictionary))
print("\n")


#save dictionary (4) state_observation_pair_dictionary. #The .npy will be added to the filename.
np.save('state_observation_pair.dict',state_observation_pair_dictionary)


"""
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$__dictionary (4)__end__$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
"""

"""
an example of dictionary item:
state_pair_dictionary['FH']
{'TTTCAT': [0, 0, 0], 'TTTCAC': [0, 0, 0], 'TTCCAT': [0, 0, 0], 'TTCCAC': [0, 0, 0]}
For the [0, 0, 0], the 1st 0 designates the feature weight, 2nd the empirical count, 3rd the expected count.

The situation of state_observation_pair_dictionary is similar to state_pair_dictionary.
'GCTA': [0, 0, 0], the 1st 0 designates the feature weight, 2nd the empirical count, 3rd the expected count.
"""


###(5)obervation_states_dictionary for forward and backward
observation_states_dictionary = amino_acid_to_codons_dictionary | first_amino_acid_to_codons_dictionary
observation_states_dictionary['START'] = ['START']
observation_states_dictionary['STOP'] = ['STOP']

print("\n")
print("Here is the obervation_states_dictionary.")
print(observation_states_dictionary)
print("The length of the obervation_states_dictionary is %s." % len(observation_states_dictionary))
print("\n")


#save dictionary (5) obervation_states_dictionary. #The .npy will be added to the filename.
np.save('observation_states.dict',observation_states_dictionary)          


