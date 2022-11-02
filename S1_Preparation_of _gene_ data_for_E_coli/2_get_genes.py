"""
The py gets genes from E. coli genome and stores in the file
The genome and gene/protein files are both downloaded from NCBI.
"""

import os

#dictionary for base complement
base_complement_dictionary = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}


#get the genome sequence from file and stored in a string
f = open('MG1655.fasta','r')
genome_file_lines = f.readlines()
f.close()

genome_sequence = ''
genome_file_lines_len = len(genome_file_lines)

for i in range(1,genome_file_lines_len):
    genome_file_line = genome_file_lines[i]
    genome_file_line = genome_file_line.strip()

    genome_sequence = genome_sequence + genome_file_line

print("The length of the genome sequence is %s." % len(genome_sequence))


#get the genes from file, first store current directory and make a new directory
current_direct = os.getcwd()
os.mkdir('gene_files')


f = open('proteins_167_161521.csv','r')
gene_file_lines = f.readlines()
f.close()


gene_file_lines_len = len(gene_file_lines)
gene_number = gene_file_lines_len - 1

print("The length of the gene information file is %s." % gene_file_lines_len)
print("The number of genes is %s." % gene_number)


#get each gene sequence and write to a txt file
gene = ''
for i in range(1,gene_file_lines_len):
    gene_file_line = gene_file_lines[i]
    gene_file_line = gene_file_line.strip()
    gene_file_line_split_to_list = gene_file_line.split(',')

    start_pos_num = int(gene_file_line_split_to_list[2])
    stop_pos_num = int(gene_file_line_split_to_list[3])
    strand_sign = gene_file_line_split_to_list[4]
    gene_id = gene_file_line_split_to_list[5]   

    if strand_sign[1] == '+':              
        gene = gene + genome_sequence[start_pos_num-1:stop_pos_num]
        os.chdir('gene_files')
        gene_file_name = gene_id+'.txt'
        gene_file = open(gene_file_name,'w')
        gene_file.write(gene)
        gene_file.close()
        os.chdir(current_direct)        
        gene = ''      

    if strand_sign[1] == '-':
        gene_complement = ''
        gene = gene + genome_sequence[start_pos_num-1:stop_pos_num]
        for j in range(len(gene)):
            base = gene[j]
            base_complement = base_complement_dictionary[base]
            gene_complement = gene_complement + base_complement
            
        gene_complement_rev = gene_complement[::-1]
        
        os.chdir('gene_files')
        gene_file_name = gene_id+'.txt'
        gene_file = open(gene_file_name,'w')
        gene_file.write(gene_complement_rev)
        gene_file.close()
        os.chdir(current_direct)

        gene = ''
        gene_complement = ''
        gene_complement_rev = ''

        

