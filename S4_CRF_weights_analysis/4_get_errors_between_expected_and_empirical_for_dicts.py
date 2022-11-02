"""
This py file is used to get the errors between the expected and empirical of edges and vertices.
"""


import numpy as np
import mpmath


#load dictionaries
state_pair_dictionary = np.load('state_pair_CRF_final.dict.npy',allow_pickle=True).item()
state_observation_pair_dictionary = np.load('state_observation_pair_CRF_final.dict.npy',allow_pickle=True).item()

error_between_expected_and_empirical_for_edge_list = []
error_between_expected_and_empirical_for_vertice_list = []

#get error_between_expected_and_empirical_for_edge_list
for observation_pair in state_pair_dictionary.keys():
    for state_pair in state_pair_dictionary[observation_pair].keys():
        error_between_expected_and_empirical_for_edge = state_pair_dictionary[observation_pair][state_pair][2]/\
                                                        state_pair_dictionary[observation_pair][state_pair][1] - 1

        error_between_expected_and_empirical_for_edge_list.append(error_between_expected_and_empirical_for_edge)

#get error_between_expected_and_empirical_for_vertice_list
for state_observation in state_observation_pair_dictionary.keys()-{'STARTSTART','STOPSTOP'}:
    error_between_expected_and_empirical_for_vertice = state_observation_pair_dictionary[state_observation][2]/\
                                                       state_observation_pair_dictionary[state_observation][1] - 1

    error_between_expected_and_empirical_for_vertice_list.append(error_between_expected_and_empirical_for_vertice)


"""
#transform to float
for i in range(len(error_between_expected_and_empirical_for_edge_list)):
    error_between_expected_and_empirical_for_edge_list[i] = float(error_between_expected_and_empirical_for_edge_list[i])
"""

#plot the data figure 
import matplotlib.pyplot as plt
import matplotlib
import numpy as np



x=[i for i in range(len(error_between_expected_and_empirical_for_edge_list))]
y=error_between_expected_and_empirical_for_edge_list
plt.scatter(x, y, s=2)
plt.xlabel("No. of edge")
plt.ylabel("error of edge")
plt.show()

x=[i for i in range(len(error_between_expected_and_empirical_for_vertice_list))]
y=error_between_expected_and_empirical_for_vertice_list
plt.scatter(x, y, s=10)
plt.xlabel("No. of vertex")
plt.ylabel("error of vertex")
plt.show()


