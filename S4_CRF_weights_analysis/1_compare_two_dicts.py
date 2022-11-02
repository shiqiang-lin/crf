"""
This py file is used to compare the empiricals of the two pairs of dictionaries.
"""


import numpy as np
import mpmath



"""
Load four dictionaries, the keys of the dictionaries are the same, however, the corresponding values may be different.
Therefore, we need to compare those values to make sure that they are identical.
"""
state_pair_dictionary_0 = np.load('state_pair.dict.npy',allow_pickle=True).item()
state_pair_dictionary_1 = np.load('state_pair_CRF_final.dict.npy',allow_pickle=True).item()
state_observation_pair_dictionary_0 = np.load('state_observation_pair.dict.npy',allow_pickle=True).item()
state_observation_pair_dictionary_1 = np.load('state_observation_pair_CRF_final.dict.npy',allow_pickle=True).item()



#compare the empiricals of state pairs in the two dicts
for observation_pair in state_pair_dictionary_0.keys():
    for state_pair in state_pair_dictionary_0[observation_pair].keys():
        if state_pair_dictionary_0[observation_pair][state_pair][1] != state_pair_dictionary_1[observation_pair][state_pair][1]:
            print("The empiricals of %s %s are different for the two dicts." % (observation_pair,state_pair))
        #else:
            #print("OK")

#compare the empiricals of state observations in the two dicts
for state_observation in state_observation_pair_dictionary_0.keys()-{'STARTSTART','STOPSTOP'}:
    if state_observation_pair_dictionary_0[state_observation][1] != state_observation_pair_dictionary_1[state_observation][1]:
        print("The empiricals of %s are different for the two dicts." % state_observation)
