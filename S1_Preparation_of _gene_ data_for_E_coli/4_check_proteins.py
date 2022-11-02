"""
The py file checks if there is any stop codon in the middle of the gene.
The sequence excluding the first codon and the last codon is [3:-3].Then walk
through the sequence with three bases as a step to check for stop codon.

Also, this script gets the longest gene, which is useful for the S algorithm
during parameter learning of the CRF model.
"""

import os
from operator import itemgetter
import numpy as np
#from build_dictionaries import observation_states_dictionary



observation_states_dictionary = np.load('observation_states.dict.npy',allow_pickle=True).item()



current_direct = os.getcwd()

gene_file_names_list = os.listdir('gene_files')
gene_file_names_list_len = len(gene_file_names_list)

collected_start_codons = []
collected_stop_codons = []

gene_id_lengths = []
gene_id_length = []


os.chdir('gene_files')
gene_file_name = ''
for i in range(gene_file_names_list_len):   
    gene_file_name = gene_file_names_list[i]


    #file name must ended with 'txt'
    gene_file_name_split = gene_file_name.strip().split('.')
    
    if gene_file_name_split[-1] == 'txt':    
        gene_file_name_f = open(gene_file_name)
        gene_file_name_content = gene_file_name_f.read()
        gene_file_name_f.close()

        gene_file_name_content = gene_file_name_content.strip()

        #record gene file name and length(bps)
        gene_id_length.append(gene_file_name)
        gene_id_length.append(len(gene_file_name_content))
        
        gene_id_lengths.append(gene_id_length)
        gene_id_length = []
            
        gene_file_name_content_middle = gene_file_name_content[3:-3]

        if len(gene_file_name_content_middle) % 3 != 0:
            print("The gene %s has a length that is not the folds of three.\n" % gene_file_name)
        
        codon_num = int(len(gene_file_name_content_middle)/3)

        for j in range(codon_num):
            codon = gene_file_name_content_middle[3*j:3*j+3]

            if codon in ['TAA', 'TGA', 'TAG']:
                print("The gene %s contains at least one internal stop codon." % gene_file_name)

                position_begin,position_end = 3*j+1+3,3*j+3+3
                """
                in position_begin:
                    +1, the base starts at 1;
                    +3, consider the start codon now;
                in posision_end:
                    +3,two bases more after the _begin position;
                    +3, consider the start codon now;
                """            
                print("The stop codon position in the gene is %s-%s:%s.\n" % (position_begin,position_end,codon))            
                break
else:
    print("\n")
    print("The proteins are all okay.\n")

    
#sort the gene_id_lengths
gene_id_lengths_sorted = sorted(gene_id_lengths, key=itemgetter(1))
print("Here are the longest gene and its length(bps)")
print(gene_id_lengths_sorted[-1])

print("Here are the shortest gene and its length(bps)")
print(gene_id_lengths_sorted[0])


      
