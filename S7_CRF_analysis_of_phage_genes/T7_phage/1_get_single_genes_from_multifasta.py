"""
This py file is used to get single genes from a multifasta file.
The single gene txt files are stored in a new folder.
"""
import os
import shutil
from Bio import SeqIO


current_direct = os.getcwd()
file_names_list = os.listdir(current_direct)

#make a dirctory gene_files, if exists, delete and make the dirctory
foldername = 'gene_files'
dirs = os.listdir(current_direct)

if(foldername not in dirs):
    os.mkdir(foldername)
else:
    shutil.rmtree(foldername)
    os.mkdir(foldername)


#print the gene sequence to a file
for seq_record in SeqIO.parse("T7_genes.txt","fasta"):
    os.chdir(foldername)
    
    gene_file_name = open(seq_record.id + '.txt','w')
    print(seq_record.seq,file=gene_file_name)
    gene_file_name.close()
    os.chdir(current_direct)



