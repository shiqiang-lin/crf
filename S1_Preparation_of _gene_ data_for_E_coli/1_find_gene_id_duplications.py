"""
The py file is used to find the duplication of gene ids.
"""
from collections import Counter


#open the protein record file and get the records
f = open('proteins_167_161521.csv','r')
gene_file_lines = f.readlines()
f.close()

gene_file_lines_len = len(gene_file_lines)
gene_number = gene_file_lines_len - 1

print("The length of the gene information file is %s." % gene_file_lines_len)
print("The number of genes is %s." % gene_number)


gene_id_list = []
for i in range(1,gene_file_lines_len):
    gene_file_line = gene_file_lines[i]
    gene_file_line = gene_file_line.strip()
    gene_file_line_split_to_list = gene_file_line.split(',')
   
    gene_id_list.append(gene_file_line_split_to_list[5])

print("The length of list:")
print(len(gene_id_list))
print("The length of set:")
print(len(set(gene_id_list)))

#get the duplicated IDs
b_dict = dict(Counter(gene_id_list))
print("Here are the duplications:")
print([key for key,value in b_dict.items() if value > 1])  #print duplicate elements
print({key:value for key,value in b_dict.items() if value > 1})  #print duplicate elements and times



