"""
This py file is used to get the codon identity percentage for each pair of the gene sequence and the max
gene sequence. The gene names and the codon identity percentages are sorted and stored in a new txt file.
"""

import os
from operator import itemgetter
import numpy as np


#load dictionary
possibility_of_identity_states_without_preference = np.load('possibility_of_identity_states_without_preference.dict.npy',allow_pickle=True).item()


current_directory = os.getcwd()

file_names_list = os.listdir('gene_files')

gene_files_list = []
#use txt in the file name to make sure it is a gene file
for file_name in file_names_list:
    file_name_split = file_name.strip().split('.')
    if file_name_split[-1] == 'txt':
        gene_files_list.append(file_name)


file_max_names_list = os.listdir('gene_max_files')

gene_max_files_list = []
#use txt in the file name to make sure it is a gene_max file
for file_name in file_max_names_list:
    file_name_split = file_name.strip().split('.')
    if file_name_split[1] == 'txt':
        gene_max_files_list.append(file_name)


"""
The input is two gene sequences and the output is ratio of identical codons
at the same positions for the two random gene sequences. 
"""
def get_ratio_of_identical_codons_from_two_gene_sequences(s0,s1):   
    number_of_identical_codons = 0
    number_of_total_codons = int(len(s1)/3)
    for i in range(number_of_total_codons):
        codon_0 = s0[3*i:3*i+3]
        codon_1 = s1[3*i:3*i+3]
        if codon_0 == codon_1:
            number_of_identical_codons = number_of_identical_codons + 1
        ratio_of_identical_codons = number_of_identical_codons/number_of_total_codons

    return ratio_of_identical_codons



ratios_of_identical_codons_of_genes_list = []
ratio_of_identical_codons_of_genes_list = []

for gene_file in gene_files_list:
    os.chdir(current_directory)
    os.chdir('gene_files')
    f = open(gene_file,'r')
    s0 = f.read().strip()
    f.close()

    gene_file_split_list = gene_file.split('.')
    gene_max_file = ''.join(gene_file_split_list[0:-1]) + '_max.txt'

    os.chdir(current_directory)
    os.chdir('gene_max_files')
    f = open(gene_max_file,'r')
    s1 = f.read().strip()
    f.close()

  

    #enter the protein files folder to get each protein file
    os.chdir(current_directory)
    #load dictionary
    counting_observations_of_protein_sequence = np.load('counting_observations_of_protein_sequence.dict.npy',allow_pickle=True).item()
    os.chdir('protein_files')
    f = open(gene_file,'r')
    protein_sequence = f.read().strip()
    f.close()

    for i in range(len(protein_sequence)):
        if i == 0 and protein_sequence[i] != 'M':
            observation = protein_sequence[i] + 's'
        else:
            observation = protein_sequence[i] 
        
        counting_observations_of_protein_sequence[observation] = counting_observations_of_protein_sequence[observation] + 1

    sum_of_counts = 0
    for observation in counting_observations_of_protein_sequence.keys():
        sum_of_counts = sum_of_counts + counting_observations_of_protein_sequence[observation]
    for observation in counting_observations_of_protein_sequence.keys():   
        counting_observations_of_protein_sequence[observation] = counting_observations_of_protein_sequence[observation]/sum_of_counts
    

    """
    To get the mathematical expectation of identity percentage, we need to sum up each observation percentage multiplied by the
    possibility of the same state for this observation.
    """
    mathematical_expectation_of_identity_percentage_without_preference = 0
    for observation in counting_observations_of_protein_sequence.keys():
        mathematical_expectation_of_identity_percentage_without_preference = \
                                                mathematical_expectation_of_identity_percentage_without_preference \
                                                + possibility_of_identity_states_without_preference[observation] \
                                                * counting_observations_of_protein_sequence[observation]

    ratio_of_identical_codons_from_two_gene_sequences = get_ratio_of_identical_codons_from_two_gene_sequences(s0,s1)
    ratio_of_identical_codons_of_genes_list.append(gene_file)
    ratio_of_identical_codons_of_genes_list.append(gene_max_file)
    ratio_of_identical_codons_of_genes_list.append(ratio_of_identical_codons_from_two_gene_sequences)
    
    ratio_of_identical_codons_of_genes_list.append(mathematical_expectation_of_identity_percentage_without_preference)
    ratio_of_identical_codons_of_genes_list.append\
                (ratio_of_identical_codons_from_two_gene_sequences-mathematical_expectation_of_identity_percentage_without_preference)

    ratios_of_identical_codons_of_genes_list.append(ratio_of_identical_codons_of_genes_list)
    ratio_of_identical_codons_of_genes_list = []
"""
[0] gene_file, [1] gene_max_file, [2] identity_percentage between [0] and [1],
[3] mathematical_expectation_of_identity_percentage_without_preference,
[4] ([2]-[3])
"""
ratios_of_identical_codons_of_genes_list_sorted = sorted(ratios_of_identical_codons_of_genes_list,key = itemgetter(4,0))


#write file
os.chdir(current_directory)

file = open('ratios_of_identical_codons_of_genes.txt','w')

for i in range(len(ratios_of_identical_codons_of_genes_list_sorted)):
    file.write(str(ratios_of_identical_codons_of_genes_list_sorted[i])+'\r\n')
file.close()







#plot the figures
import matplotlib.pyplot as plt


difference_between_gene_gene_max_list = []
mathematical_expectation_of_identity_percentage_without_preference_list = []

for i in range(len(ratios_of_identical_codons_of_genes_list_sorted)):
    difference_between_gene_gene_max_list.append(ratios_of_identical_codons_of_genes_list_sorted[i][2])
    mathematical_expectation_of_identity_percentage_without_preference_list.append(ratios_of_identical_codons_of_genes_list_sorted[i][3])
x = [i for i in range(len(difference_between_gene_gene_max_list))]
y1 = difference_between_gene_gene_max_list
y2 = mathematical_expectation_of_identity_percentage_without_preference_list

plt.scatter(x,y1,c='royalblue',s=20)
plt.scatter(x,y2,c='darkorange',s=20)

plt.xlabel("gene")
plt.ylabel("codon usage index")
plt.title("phage T4")
plt.show()


difference_between_gene_gene_max_and_mathematical_expectation_without_preference_list = []
for i in range(len(ratios_of_identical_codons_of_genes_list_sorted)):
    difference_between_gene_gene_max_and_mathematical_expectation_without_preference_list.append\
                                                    (ratios_of_identical_codons_of_genes_list_sorted[i][4])
    
plt.hist(difference_between_gene_gene_max_and_mathematical_expectation_without_preference_list,10)
plt.xlabel("difference of codon usage index")
plt.ylabel("frequency number")
plt.title("phage T4")
plt.show()





    
