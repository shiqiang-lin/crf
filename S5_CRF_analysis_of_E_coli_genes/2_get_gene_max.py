"""
This py file is used to get the gene sequence of max probability for each protein sequence of E. coli
with Viterbi algorithm. The max gene sequences are stored in a new folder.
"""

import os
import shutil
import numpy as np
import mpmath
from operator import itemgetter



current_directory = os.getcwd()

#load dictionaries
observation_states_dictionary = np.load('observation_states.dict.npy',allow_pickle=True).item()
state_pair_dictionary = np.load('state_pair_CRF_final.dict.npy',allow_pickle=True).item()
state_observation_pair_dictionary = np.load('state_observation_pair_CRF_final.dict.npy',allow_pickle=True).item()


#make a dirctory gene_max_files, if exists, delete and make the dirctory
foldername = 'gene_max_files'
dirs = os.listdir(current_directory)

if(foldername not in dirs):
    os.mkdir(foldername)
else:
    shutil.rmtree(foldername)
    os.mkdir(foldername)

file_names_list = os.listdir('gene_files')

protein_files_list = []
#use txt in the file name to make sure it is a protein file
for file_name in file_names_list:
    file_name_split = file_name.strip().split('.')
    if file_name_split[1] == 'txt':
        protein_files_list.append(file_name)


"""
The input is a string of protein sequence and the output is a string of DNA sequence.
"""
def get_max_gene_sequence(protein_sequence):
    #get observation_list and first add 'START'
    observation_list = []
    observation_list.append('START')

    for i in range(len(protein_sequence)):
        if i == 0 and protein_sequence[i] != 'M':
            observation = protein_sequence[i] + 's'
        else:
            observation = protein_sequence[i]
        
        observation_list.append(observation)
    observation_list.append('STOP')

    #initialize Viterbi
    right_state_delta_record_phi_record_list = [{'START':{'delta_record':0, 'phi_record':'START'}}]

    #set tmp containers
    left_state_right_state_exponent_list = []
    left_state_right_state_exponents_list = []
    delta_record_phi_record_dictionary = {}
    right_state_delta_record_phi_record_dictionary = {}
    
    
    #forward
    for i in range(0,len(observation_list)-1):
        left_observation = observation_list[i]        
        left_states = observation_states_dictionary[left_observation]
        right_observation = observation_list[i+1]        
        right_states = observation_states_dictionary[right_observation]
        
        for right_state in right_states:
            for left_state in left_states:
                left_state_right_state_exponent_list.append(left_state)                    
                left_state_right_state_exponent_list.append \
                                    (right_state_delta_record_phi_record_list[-1][left_state]['delta_record'] + \
                                    float(state_pair_dictionary[left_observation+right_observation][left_state+right_state][0]) + \
                                    float(state_observation_pair_dictionary[right_state+right_observation][0]))
                left_state_right_state_exponents_list.append(left_state_right_state_exponent_list)
                left_state_right_state_exponent_list = []
            left_state_right_state_exponents_list_sorted = sorted(left_state_right_state_exponents_list,key=itemgetter(1))

            delta_record_phi_record_dictionary['delta_record'] = left_state_right_state_exponents_list_sorted[-1][1]
            delta_record_phi_record_dictionary['phi_record'] = left_state_right_state_exponents_list_sorted[-1][0]

            left_state_right_state_exponents_list = []
            left_state_right_state_exponents_list_sorted = []
            

            right_state_delta_record_phi_record_dictionary[right_state] = delta_record_phi_record_dictionary
            delta_record_phi_record_dictionary = {}


        right_state_delta_record_phi_record_list.append(right_state_delta_record_phi_record_dictionary)
        

        right_state_delta_record_phi_record_dictionary = {}


    #get the max gene sequence
    max_gene_sequence_list = []
    previous_state = 'STOP'        
    previous_previous_state = right_state_delta_record_phi_record_list[-1][previous_state]['phi_record']

    max_gene_sequence_list.append(previous_previous_state)

    for i in range(len(right_state_delta_record_phi_record_list)-2,1,-1):
        previous_state = previous_previous_state
        previous_previous_state = right_state_delta_record_phi_record_list[i][previous_state]['phi_record']
        max_gene_sequence_list.append(previous_previous_state)


    max_gene_sequence = ''
    for i in range(len(max_gene_sequence_list)-1,-1,-1):
        max_gene_sequence = max_gene_sequence + max_gene_sequence_list[i]

    return(max_gene_sequence)
        
            
#go through each protein file in the protein_files_list
for protein_file in protein_files_list:
    #open the protein sequence file
    os.chdir(current_directory)
    os.chdir('protein_files')
    
    protein_sequence_file = open(protein_file,'r')
    protein_sequence = protein_sequence_file.read().strip()
    protein_sequence_file.close()


    max_gene_sequence = get_max_gene_sequence(protein_sequence)

    os.chdir(current_directory)
    os.chdir(foldername)


    protein_file_split_list = protein_file.split('.')
    #notice that the protein file name should be named by gene ID
    g = open(protein_file_split_list[0]+'_max.'+protein_file_split_list[1],'w')
    g.write(max_gene_sequence)
    g.close()        
            
                                     

