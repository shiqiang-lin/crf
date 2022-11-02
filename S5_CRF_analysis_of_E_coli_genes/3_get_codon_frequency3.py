"""
This py file is used to get the codon frequency from the state_observation_pair_dictionary. The 
ecoli_observation_codon_frequency_dictionary_for_choices_use is saved as a dictionary file for
later use.
"""

import numpy as np

#load state_observation_pair_dictionary
state_observation_pair_dictionary = np.load('state_observation_pair_CRF_final.dict.npy',allow_pickle=True).item()
observation_states_dictionary = np.load('observation_states.dict.npy',allow_pickle=True).item()


ecoli_observation_codon_frequency_dictionary = {}
codon_frequency_dictionary = {}
for observation in observation_states_dictionary.keys()-{'START', 'STOP'}:
    for state in observation_states_dictionary[observation]:
        codon_frequency_dictionary[state] = 0
    ecoli_observation_codon_frequency_dictionary[observation] = codon_frequency_dictionary
    codon_frequency_dictionary = {}

#print(ecoli_observation_codon_frequency_dictionary)

#get the count of each observation state
for observation in ecoli_observation_codon_frequency_dictionary.keys():
    for state in ecoli_observation_codon_frequency_dictionary[observation]:
        ecoli_observation_codon_frequency_dictionary[observation][state] = state_observation_pair_dictionary[state+observation][1]

#print(ecoli_observation_codon_frequency_dictionary)
    
#get the frequency of each observation state
sum_of_counts = 0
for observation in ecoli_observation_codon_frequency_dictionary.keys():
    for state in ecoli_observation_codon_frequency_dictionary[observation]:
        sum_of_counts = sum_of_counts + ecoli_observation_codon_frequency_dictionary[observation][state]
    for state in ecoli_observation_codon_frequency_dictionary[observation]:
        ecoli_observation_codon_frequency_dictionary[observation][state] = \
                                                ecoli_observation_codon_frequency_dictionary[observation][state]/sum_of_counts
    sum_of_counts = 0
    
#print(ecoli_observation_codon_frequency_dictionary)



ecoli_observation_codon_frequency_dictionary_for_choices_use = {}
codons_list = []
frequencies_list = []
codons_frequencies_list = []
for observation in ecoli_observation_codon_frequency_dictionary.keys():
    for state in ecoli_observation_codon_frequency_dictionary[observation]:
        codons_list.append(state)
        frequencies_list.append(ecoli_observation_codon_frequency_dictionary[observation][state])
    codons_frequencies_list.append(codons_list)
    codons_frequencies_list.append(frequencies_list)

    ecoli_observation_codon_frequency_dictionary_for_choices_use[observation] = codons_frequencies_list
    codons_list = []
    frequencies_list = []
    codons_frequencies_list = []

print(ecoli_observation_codon_frequency_dictionary_for_choices_use)


np.save('ecoli_observation_codon_frequency_dictionary_for_choices_use.dict',\
        ecoli_observation_codon_frequency_dictionary_for_choices_use)    

"""
dictionary of the possibility of identical states for observation, used to get the mathematical expectation of
percentage of the identical codons in the protein
If for one particular observation, the possibility of state_1 is p1,the possibility of state_2 is p2, ... , the
possibility of state_k is pk, then the possibility of the same state for this observation is p1*p1 + p2*p2
+ ... + pk*pk.
"""
possibility_of_identity_states = {}
for observation in ecoli_observation_codon_frequency_dictionary_for_choices_use.keys():
    sum_of_squares_of_state_possibilities = 0
    for state_possibility in ecoli_observation_codon_frequency_dictionary_for_choices_use[observation][1]:
        sum_of_squares_of_state_possibilities = sum_of_squares_of_state_possibilities + state_possibility * state_possibility
    possibility_of_identity_states[observation] = sum_of_squares_of_state_possibilities
    
print(possibility_of_identity_states)

np.save('possibility_of_identity_states.dict', possibility_of_identity_states)


possibility_of_identity_states_without_preference = {}
for observation in ecoli_observation_codon_frequency_dictionary_for_choices_use.keys():
    euqal_possibility = 1/len(ecoli_observation_codon_frequency_dictionary_for_choices_use[observation][0])   
    possibility_of_identity_states_without_preference[observation] = euqal_possibility

print(possibility_of_identity_states_without_preference)
np.save('possibility_of_identity_states_without_preference.dict', possibility_of_identity_states_without_preference)
    


#dictionary for counting observations of each protein sequence, needed when calculating mathematic expectation
counting_observations_of_protein_sequence = {}
for observation in possibility_of_identity_states.keys():
    counting_observations_of_protein_sequence[observation] = 0

np.save('counting_observations_of_protein_sequence.dict', counting_observations_of_protein_sequence)



