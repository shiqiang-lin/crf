"""
This py file is used to get the numbers of edges and vertices.
"""


import numpy as np

#load dictionaries
state_pair_dictionary = np.load('state_pair_CRF_final.dict.npy',allow_pickle=True).item()
state_observation_pair_dictionary = np.load('state_observation_pair_CRF_final.dict.npy',allow_pickle=True).item()

number_of_edges = 0
number_of_vertices = 0

for observation_pair in state_pair_dictionary.keys():
    number_of_edges = number_of_edges + len(state_pair_dictionary[observation_pair])

print("The number of edges is %s." % number_of_edges)

number_of_vertices = number_of_vertices + len(state_observation_pair_dictionary)

print("The number of vertices is %s." % number_of_vertices)


