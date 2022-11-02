"""
This py file is used to get each protein sequence from the corresponding gene sequence
file in the folder.
"""

import os
import shutil

#translate the DNA sequence to protein sequence
def translate_gene_to_protein(gene_sequence):
    codon_dict = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TGT': 'C',\
                    'TGC': 'C', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',\
                    'TTT': 'F', 'TTC': 'F', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',\
                    'GGG': 'G', 'CAT': 'H', 'CAC': 'H', 'ATT': 'I', 'ATC': 'I',\
                    'ATA': 'I', 'AAA': 'K', 'AAG': 'K', 'TTA': 'L', 'TTG': 'L',\
                    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATG': 'M',\
                    'AAT': 'N', 'AAC': 'N', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',\
                    'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',\
                    'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'TCT': 'S',\
                    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',\
                    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GTT': 'V',\
                    'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TGG': 'W', 'TAT': 'Y',\
                    'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGA': '*'
                  }

    protein = ''
    for i in range(int(len(gene_sequence)/3)):
        protein = protein + codon_dict[gene_sequence[3*i:3*i+3]]

    return protein
    """
    print("\n",end='')       
    print("The translated amino acid sequence of the protein is ")  
    print(protein)
    """

current_direct = os.getcwd()
file_names_list = os.listdir('gene_files')

#make a dirctory protein_files, if exists, delete and make the dirctory
foldername = 'protein_files'
dirs = os.listdir(current_direct)

if(foldername not in dirs):
    os.mkdir(foldername)
else:
    shutil.rmtree(foldername)
    os.mkdir(foldername)


gene_files_list = []
#use txt in the file name to make sure it is a gene file
for file_name in file_names_list:
    file_name_split = file_name.strip().split('.')
    if file_name_split[1] == 'txt':
        gene_files_list.append(file_name)
        
#get protein sequence and store in txt file 
for i in range(len(gene_files_list)):
    gene_file = gene_files_list[i]
    os.chdir('gene_files')
    f = open(gene_file,'r')
    gene_sequence = f.read().strip()
    f.close()
    
    translated_protein_sequence = translate_gene_to_protein(gene_sequence)
    os.chdir(current_direct)
    os.chdir('protein_files')

    #notice that the protein file name should be named by gene ID
    #and should not use the gene_file (line 53) above. 
    g = open(gene_file,'w')
    g.write(translated_protein_sequence)
    g.close()
    os.chdir(current_direct)


