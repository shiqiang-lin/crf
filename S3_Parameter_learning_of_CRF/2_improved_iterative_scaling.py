"""
This py file conducts forward algorithm, backward algorithm, and feature expectatoin calculation for each gene.
The feature counts in the two dictionaries are updated and the weight of each feature is updated with the slack
feature. The S algorithm of the Improved Iterative Scaling (IIS)is used to optimize the weight of each feature
in the two dictionaries.
"""

import os
import sys
import datetime
from multiprocessing import Process,set_start_method
import numpy as np
from mpmath import exp
from mpmath import log



current_directory = os.getcwd()
#os.chdir(current_directory)

#load dictionaries
codon_to_amino_acid_dictionary = np.load('codon_to_amino_acid.dict.npy',allow_pickle=True).item()
amino_acid_to_codons_dictionary = np.load('amino_acid_to_codons.dict.npy',allow_pickle=True).item()
observation_states_dictionary = np.load('observation_states.dict.npy',allow_pickle=True).item()


global state_pair_dictionary
global state_observation_pair_dictionary

state_pair_dictionary = np.load('state_pair.dict.npy',allow_pickle=True).item()
state_observation_pair_dictionary = np.load('state_observation_pair.dict.npy',allow_pickle=True).item()



def prepare_genome_observation_state_data():
    os.chdir(current_directory) 
    

    #list for storing genome data of observations and states
    genome_observation_state_list = []
    gene_observation_state_list = []
    
    file_names_list = os.listdir('gene_files')

    #get only the .txt files to list
    gene_file_names_list = []   
    for file_name in file_names_list:
        file_name_split = file_name.split('.')
        if file_name_split[1] == 'txt':
            gene_file_names_list.append(file_name)

    gene_file_names_list_len = len(gene_file_names_list)

    
    os.chdir('gene_files')
    for i in range(gene_file_names_list_len):
    #for i in range(10):
        gene_file_name = gene_file_names_list[i]
        gene_file_name_f = open(gene_file_name)
        gene_file_name_content = gene_file_name_f.read()
        gene_file_name_f.close()

        gene_file_name_content = gene_file_name_content.strip()

        
        #get observation_list and first add 'START'
        observation_list = []
        observation_list.append('START')
               
        """
        Then deal with the middle part of the gene, which consists of codon_1, codon_2, ..., codon_n.

        codon_1      0,1,2
        codon_2      3,4,5
        codon_3      6,7,8
        codon_4      9,10,11
        codon_5      12,13,14
        ...          ...
        codon_i      3*i-3,3*i-2,3*i-1
                     [3*i-3,3*i]
        ...          ...
        codon_n-1    3*(n-1)-3,3*(n-1)-3+1,3*(n-1)-3+2
                     [3*n-6,3*n-3]
        codon_n      3*n-3,3*n-2,3*n-1
                     [3*n-3,3*n]          
        """                                                         
        gene_file_name_content_len = len(gene_file_name_content)
        codon_num = int(gene_file_name_content_len/3)

        for i in range(1,codon_num+1):
            codon = gene_file_name_content[3*i-3:3*i]
        #the first codon is different from the rest
            if (i == 1) and (codon != 'ATG'):
                observation = codon_to_amino_acid_dictionary[codon] + 's'
            else:
                observation = codon_to_amino_acid_dictionary[codon]                
            observation_list.append(observation)
            
        #finally add 'STOP'    
        observation_list.append('STOP')
        

        
        #lists for storing forward result and backward result
        forward_observation_states_list = []
        backward_observation_states_list = []

        #deal with forward_observation_states_list
        position_states_dictionary = {}
        position_observation_states_dictionary = {}
        
        for i in range(len(observation_list)):
            observation = observation_list[i]
            states = observation_states_dictionary[observation]
            
            for state in states:
                position_states_dictionary[state] = 0

            position_observation_states_dictionary[observation] = position_states_dictionary
            forward_observation_states_list.append(position_observation_states_dictionary)            

            #empty the containers
            position_states_dictionary = {}
            position_observation_states_dictionary = {}
            
        #deal with backward_observation_states_list
        for i in range(len(observation_list)):
            observation = observation_list[i]
            states = observation_states_dictionary[observation]
            
            for state in states:
                position_states_dictionary[state] = 0

            position_observation_states_dictionary[observation] = position_states_dictionary
            backward_observation_states_list.append(position_observation_states_dictionary)

            #empty the containers
            position_states_dictionary = {}
            position_observation_states_dictionary = {}


        #store each gene data to genome list
        gene_observation_state_list.append(forward_observation_states_list)
        gene_observation_state_list.append(backward_observation_states_list)
        genome_observation_state_list.append(gene_observation_state_list)
        #empty the containers
        gene_observation_state_list = []

    return genome_observation_state_list



def get_expected_counts(task_id):
    sub_task_id = task_id
       
    for gene_data in genome_data[task_id::sub_process_count]:
        forward_observation_states_list = gene_data[0]
        backward_observation_states_list = gene_data[1]
               
            
        #forward
        #set the state of 'START' to 1       
        forward_observation_states_list[0]['START']['START'] = 1  
        for i in range(0,len(forward_observation_states_list)-1):
            left_observation = [item for item in forward_observation_states_list[i]][0]        
            left_states = set(forward_observation_states_list[i][left_observation].keys())
            right_observation = [item for item in forward_observation_states_list[i+1]][0]
            right_states = set(forward_observation_states_list[i+1][right_observation].keys())
            
            for right_state in right_states:
                for left_state in left_states:
                    forward_observation_states_list[i+1][right_observation][right_state] = \
                        forward_observation_states_list[i+1][right_observation][right_state] + \
                        forward_observation_states_list[i][left_observation][left_state] * \
                        exp(state_pair_dictionary[left_observation+right_observation][left_state+right_state][0] + \
                            state_observation_pair_dictionary[right_state+right_observation][0])                           
                      
  
        #backward
        #set the state of 'STOP' to 1       
        backward_observation_states_list[-1]['STOP']['STOP'] = 1      
        for i in range(len(backward_observation_states_list)-1,0,-1):            
            right_observation = [item for item in backward_observation_states_list[i]][0]
            right_states = set(backward_observation_states_list[i][right_observation].keys())
            left_observation = [item for item in backward_observation_states_list[i-1]][0]
            left_states = set(backward_observation_states_list[i-1][left_observation].keys())
                     
            for left_state in left_states:
                for right_state in right_states:                    
                    backward_observation_states_list[i-1][left_observation][left_state] = \
                        backward_observation_states_list[i-1][left_observation][left_state] + \
                        backward_observation_states_list[i][right_observation][right_state] * \
                        exp(state_pair_dictionary[left_observation+right_observation][left_state+right_state][0] + \
                            state_observation_pair_dictionary[right_state+right_observation][0])                                        
      

        #caculate expected count
        #average the two sums to get nomalization. Acturally, the two sums are almost the same, with a very small difference.
        nomalization = (forward_observation_states_list[-1]['STOP']['STOP'] + \
                        backward_observation_states_list[0]['START']['START']) / 2
                        
        for i in range(0,len(forward_observation_states_list)-1):
            left_observation = [item for item in forward_observation_states_list[i]][0]        
            left_states = set(forward_observation_states_list[i][left_observation].keys())
            right_observation = [item for item in forward_observation_states_list[i+1]][0]
            right_states = set(forward_observation_states_list[i+1][right_observation].keys())

            #add state_pair_expectation to state_observation_pair_dictionary
            for left_state in left_states:
                for right_state in right_states:
                    state_pair_exponent = state_pair_dictionary[left_observation+right_observation][left_state + right_state][0] + \
                                          state_observation_pair_dictionary[right_state+right_observation][0]
                    state_pair_mass = exp(state_pair_exponent) 
                    left_mass_forward = forward_observation_states_list[i][left_observation][left_state]
                    right_mass_backward = backward_observation_states_list[i+1][right_observation][right_state]
                                      
                    state_pair_expectation = left_mass_forward * state_pair_mass * right_mass_backward / nomalization
                    state_pair_dictionary[left_observation+right_observation][left_state + right_state][2] = \
                            state_pair_dictionary[left_observation+right_observation][left_state + right_state][2] + \
                            state_pair_expectation

            #add state_observation_expectation to state_observation_pair_dictionary                 
            for right_state in right_states:                 
                    right_mass_forward = forward_observation_states_list[i+1][right_observation][right_state]
                    state_observation_expectation = right_mass_forward * right_mass_backward / nomalization
                    state_observation_pair_dictionary[right_state+right_observation][2] = \
                            state_observation_pair_dictionary[right_state+right_observation][2] + \
                            state_observation_expectation

    #save the dictionaries
    os.chdir(current_directory)
    
    state_pair_dictionary_file_name = 'state_pair'+str(sub_task_id)+'.dict'
    #print(state_pair_dictionary_file_name)
    state_observation_pair_dictionary_file_name = 'state_observation_pair'+str(sub_task_id)+'.dict'
    
    np.save(state_pair_dictionary_file_name,state_pair_dictionary)                         
    np.save(state_observation_pair_dictionary_file_name,state_observation_pair_dictionary)
    
                   
                                

def merge_and_revise_dictionaries():
    os.chdir(current_directory)

    #load dictionaries stored by subprocesses
    loaded_state_pair_dictionaries_list = []
    loaded_state_observation_pair_dictionaries_list = []
    for i in range(sub_process_count):
              
        state_pair_dictionary_name = np.load('state_pair'+str(i)+'.dict.npy',allow_pickle=True).item()
        loaded_state_pair_dictionaries_list.append(state_pair_dictionary_name)      
          
        state_observation_pair_dictionary_name = np.load('state_observation_pair'+str(i)+'.dict.npy',allow_pickle=True).item()
        loaded_state_observation_pair_dictionaries_list.append(state_observation_pair_dictionary_name)


    #merge values of the dictionaries stored by subprocesses
    for observation_pair in state_pair_dictionary.keys():
        for state_pair in state_pair_dictionary[observation_pair]:
            for dictionary in loaded_state_pair_dictionaries_list:
                state_pair_dictionary[observation_pair][state_pair][2] = state_pair_dictionary[observation_pair][state_pair][2] + \
                                                                      dictionary[observation_pair][state_pair][2] 
    
    for state_observation in state_observation_pair_dictionary.keys()-{'STARTSTART', 'STOPSTOP'}:
        for dictionary in loaded_state_observation_pair_dictionaries_list:
            state_observation_pair_dictionary[state_observation][2] = state_observation_pair_dictionary[state_observation][2] + \
                                                               dictionary[state_observation][2]
   

    #see if the value is zero, if is, set a small value, 1e-8
    for observation_pair,state_pair_dict in state_pair_dictionary.items():
        for state_pair,value in state_pair_dict.items():
            if value[2] == 0:       
                #print("The %s %s has an expected feature of zero." % (observation_pair,state_pair))
                state_pair_dictionary[observation_pair][state_pair][2] = 1e-8
            
           
    #This dictionary only has START, STOP with zero values, which do not need to set a small value.
    for state_observation in state_observation_pair_dictionary.keys()-{'STARTSTART','STOPSTOP'}:
        if state_observation_pair_dictionary[state_observation][2] == 0:
            state_observation_pair_dictionary[state_observation][2] = 1e-8
            #print("The %s has an expected feature of zero." % state_observation)



"""
Here are the longest gene and its length(bps)
['946498.txt', 7077]
slack_feature = count of edge + count of vetices + 1 = ((7077/3+2)-1) + (7077/3+2) + 1  
"""
delta_weight_list = []
slack_feature = int(7077/3+2)*2
cycle = 0

#function for updating the weights in the two dictionaries
def update_dictionaries():
  
    global delta_weight_list
    
    for observation_pair in state_pair_dictionary.keys():
        for state_pair in state_pair_dictionary[observation_pair]:
            delta_weight_of_state_pair = 1/slack_feature * log(state_pair_dictionary[observation_pair][state_pair][1] / \
                                         state_pair_dictionary[observation_pair][state_pair][2])
            state_pair_dictionary[observation_pair][state_pair][0] = state_pair_dictionary[observation_pair][state_pair][0] + \
                                                         delta_weight_of_state_pair
            delta_weight_list.append(abs(delta_weight_of_state_pair))            
        
    for state_observation in state_observation_pair_dictionary.keys()-{'STARTSTART', 'STOPSTOP'}:
        delta_weight_state_observation = 1/slack_feature * log(state_observation_pair_dictionary[state_observation][1] / \
                                         state_observation_pair_dictionary[state_observation][2])
        state_observation_pair_dictionary[state_observation][0] = state_observation_pair_dictionary[state_observation][0] + \
                                                                  delta_weight_state_observation
        delta_weight_list.append(abs(delta_weight_state_observation))





#get genome data
genome_data = prepare_genome_observation_state_data()

#Important!!! This parameter depends on the processes of the CPU and needs to be modified for increasing computing speed. 
sub_process_count = 72
global sub_process_list
sub_process_list = []



#set expected counts to zero
for observation_pair in state_pair_dictionary.keys():
    for state_pair in state_pair_dictionary[observation_pair]:
        state_pair_dictionary[observation_pair][state_pair][2] = 0           
for state_observation in state_observation_pair_dictionary.keys():
    state_observation_pair_dictionary[state_observation][2] = 0

#start subprocesses
set_start_method('fork')

for task_id in range(sub_process_count):
    p = Process(target=get_expected_counts,args=(task_id,))
    p.start()
    sub_process_list.append((task_id,p))
    
#wait and exit subprocesses
for (task_id, p) in sub_process_list:
    p.join()

merge_and_revise_dictionaries()
update_dictionaries()

#delete the intermedian dictionaries
for i in range(sub_process_count):  
    os.system('rm state_pair'+str(i)+'.dict.npy')
    os.system('rm state_observation_pair'+str(i)+'.dict.npy')


max_delta_weight = max(delta_weight_list)
sum_of_squares_of_delta_weight = 0


if max_delta_weight <= 1e-6:
    print("The max_delta_weight has been already less than 1e-6.")
    print(datetime.datetime.now())    
    print("The cycle number: %s" % cycle)
    print("The max_delta_weight: %s" % max_delta_weight)
    print("\n",end='')

    os.chdir(current_directory)
    #save dictionaries for later use, the dictionary name add _CRF
    np.save('state_pair_CRF_final.dict',state_pair_dictionary)                            #filename with .npy
    np.save('state_observation_pair_CRF_final.dict',state_observation_pair_dictionary)    #filename with .npy
    #show on the screen for visual comparing empircal with expected
    print(state_pair_dictionary)
    print("\n")
    print(state_observation_pair_dictionary)
        
    sys.exit()
else:
    #improved iterative scaling
    while max_delta_weight > 1e-6:

        #set expected counts to zero
        for observation_pair in state_pair_dictionary.keys():
            for state_pair in state_pair_dictionary[observation_pair]:
                state_pair_dictionary[observation_pair][state_pair][2] = 0           
        for state_observation in state_observation_pair_dictionary.keys():
            state_observation_pair_dictionary[state_observation][2] = 0


        sub_process_list = []   
        for task_id in range(sub_process_count):
            p = Process(target=get_expected_counts,args=(task_id,))
            p.start()
            sub_process_list.append((task_id,p))          
        for (task_id, p) in sub_process_list:
            p.join()
        sub_process_list = []


        delta_weight_list = []
        merge_and_revise_dictionaries()
        update_dictionaries()
        
        #delete the intermedian dictionaries
        for i in range(sub_process_count):  
            os.system('rm state_pair'+str(i)+'.dict.npy')
            os.system('rm state_observation_pair'+str(i)+'.dict.npy')

        max_delta_weight = max(delta_weight_list)

        sum_of_squares_of_delta_weight = 0
        for i in range(len(delta_weight_list)):
            sum_of_squares_of_delta_weight = sum_of_squares_of_delta_weight + delta_weight_list[i] * delta_weight_list[i]       

        cycle = cycle + 1
        print(datetime.datetime.now())    
        print("The cycle number: %s" % cycle)
        print("The max_delta_weight: %s" % max_delta_weight)
        #print("\n",end='')

        print("The sum_of_squares_of_delta_weight: %s" % sum_of_squares_of_delta_weight)
        print("\n",end='')
            
        os.chdir(current_directory)
        #save dictionaries for later use, the dictionary name add _CRF
        np.save('state_pair_running.dict',state_pair_dictionary)                            
        np.save('state_observation_pair_running.dict',state_observation_pair_dictionary)    

    os.chdir(current_directory)

    #save dictionaries for later use, the dictionary name add _CRF
    np.save('state_pair_CRF_final.dict',state_pair_dictionary)                            
    np.save('state_observation_pair_CRF_final.dict',state_observation_pair_dictionary)    
    #show on the screen for visual comparing empircal with expected
    print(state_pair_dictionary)
    print("\n")
    print(state_observation_pair_dictionary)



