"""
This py file is used to print the weights of the two dictionaries.
The output state_pair_dictionary had better be copied to TextEdit for reading.
Right click the "Sqeezed text (2213 lines)", click "copy", open the TextEdit and paste.
"""


import numpy as np



#load dictionaries
state_pair_dictionary = np.load('state_pair_CRF_final.dict.npy',allow_pickle=True).item()
state_observation_pair_dictionary = np.load('state_observation_pair_CRF_final.dict.npy',allow_pickle=True).item()


print("Here is the state_pair_dictionary.")
print(state_pair_dictionary)
print("\n")

print("Here is the state_observation_pair_dictionary.")
print(state_observation_pair_dictionary)
print("\n")
