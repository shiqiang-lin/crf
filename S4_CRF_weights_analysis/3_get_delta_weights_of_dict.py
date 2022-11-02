"""
This py file is used to get the delta weights and the weights of the two dictionaries.
"""


import numpy as np
import mpmath
from mpmath import log


#load dictionaries
state_pair_dictionary = np.load('state_pair_CRF_final.dict.npy',allow_pickle=True).item()
state_observation_pair_dictionary = np.load('state_observation_pair_CRF_final.dict.npy',allow_pickle=True).item()

delta_weight_of_edge_list = []
delta_weight_of_vertice_list = []
weight_of_edge_list = []
weight_of_vertice_list = []

slack_feature = 4722
#get delta_weight_of_edge_list, weight_of_edge_list = []
for observation_pair in state_pair_dictionary.keys():
    for state_pair in state_pair_dictionary[observation_pair].keys():
        delta_weight_of_edge = log(state_pair_dictionary[observation_pair][state_pair][1]/\
                                                        state_pair_dictionary[observation_pair][state_pair][2])/4722
        delta_weight_of_edge_list.append(delta_weight_of_edge)
        weight_of_edge_list.append(state_pair_dictionary[observation_pair][state_pair][0])
        #look for values less than -10.0 and print
        if state_pair_dictionary[observation_pair][state_pair][0] < -10.0:
            print(observation_pair,state_pair,state_pair_dictionary[observation_pair][state_pair])

#get delta_weight_of_vertice_list, weight_of_vertice_list
for state_observation in state_observation_pair_dictionary.keys()-{'STARTSTART','STOPSTOP'}:
    delta_weight_of_vertice = log(state_observation_pair_dictionary[state_observation][1]/\
                                                       state_observation_pair_dictionary[state_observation][2])/4722
    delta_weight_of_vertice_list.append(delta_weight_of_vertice)
    weight_of_vertice_list.append(state_observation_pair_dictionary[state_observation][0])
    #look for values less than -2.0 and print
    if state_observation_pair_dictionary[state_observation][0] < -1.0:
        print(state_observation,state_observation_pair_dictionary[state_observation])
   


"""
#transform to float
for i in range(len(error_between_expected_and_empirical_for_edge_list)):
    error_between_expected_and_empirical_for_edge_list[i] = float(error_between_expected_and_empirical_for_edge_list[i])
"""

#plot the four data figures 
import matplotlib.pyplot as plt
import matplotlib
import numpy as np



x=[i for i in range(len(delta_weight_of_edge_list))]
y=delta_weight_of_edge_list
plt.scatter(x,y,s=2)
plt.xlabel("No. of edge")
plt.ylabel("delta weight of edge")
plt.show()


x=[i for i in range(len(delta_weight_of_vertice_list))]
y=delta_weight_of_vertice_list
plt.scatter(x,y,s=10)
plt.xlabel("No. of vertex")
plt.ylabel("delta weight of vertex")
plt.show()


#x=[i for i in range(len(weight_of_edge_list))]
#y=weight_of_edge_list
#transform to float
for i in range(len(weight_of_edge_list)):
    weight_of_edge_list[i] = float(weight_of_edge_list[i])

plt.hist(weight_of_edge_list,10)
plt.xlabel("weight of edge")
plt.ylabel("frequency number")
plt.show()


#x=[i for i in range(len(weight_of_vertice_list))]
#y= weight_of_vertice_list
#transform to float
for i in range(len(weight_of_vertice_list)):
    weight_of_vertice_list[i] = float(weight_of_vertice_list[i])

plt.hist(weight_of_vertice_list,10)
plt.xlabel("weight of vertex")
plt.ylabel("frequency number")
plt.show()




