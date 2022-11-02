"""
The py file is used to check all the genes extracted from the genome.
(1) check if all the characters are 'ATCG';
(2) get all the start codons and compare with those in the state pair dictionary;
(3) get all the stop codons and compare with those in the state pair dictionary;
(4) check if the number of characters are folds of three;
"""


import os
from collections import Counter



current_direct = os.getcwd()

gene_file_names_list = os.listdir('gene_files')
gene_file_names_list_len = len(gene_file_names_list)

collected_start_codons = []
collected_stop_codons = []


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
        gene_file_name_content_len = len(gene_file_name_content)


        collected_start_codons.append(gene_file_name_content[0:3])
        """
        delete genes (only 3) containg start codons that the E. coli does not have. This makes script works while keeps
        the result almost the same.
        """
        
        if gene_file_name_content[0:3] == 'ATA' or gene_file_name_content[0:3] == 'ATC':
            print(gene_file_name)



        
        collected_stop_codons.append(gene_file_name_content[-3:])

        if set(gene_file_name_content) != {'A', 'T', 'C', 'G'}:
            print("The gene %s contains characters other than ATCG." % gene_file_name)

        if gene_file_name_content_len % 3 != 0:
            print("The length of the gene %s is not the folds of three." % gene_file_name)
os.chdir(current_direct)
        

#get start codons and counts, stop codons and counts
start_codon_dict = dict(Counter(collected_start_codons))
print("Here are the start codons and numbers:")
for key in start_codon_dict.keys():
    print(key,end=' ')
print('\n',end='')
print({key:value for key,value in start_codon_dict.items()},'\n')


stop_codon_dict = dict(Counter(collected_stop_codons))
print("Here are the stop codons and numbers:")
for key in stop_codon_dict.keys():
    print(key,end=' ')
print('\n',end='')
print({key:value for key,value in stop_codon_dict.items()})




