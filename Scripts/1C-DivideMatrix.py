#!/usr/bin/env python3
#This script divides an expression matrix file into organ files, using an annotation file as reference
#Erielle Villanueva - 30 Aug 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/Organ/'

#Function to create a file, function takes in name of organ and set of experiments as arguments
def file_creator(org, exps):
  my_list = ['gene\t' + '\t'.join(exps) + '\n'] #creates first row of the file
  for i in range(len(gene_col)):
    row = gene_col[i] #first word is gene
    for exp in exps:
      row += '\t' + exp_dict[exp][i] #adds tpm values of exps
    row += '\n'
    my_list.append(row)
  v = open(path + 'Organ/%s.txt' % org,'w')
  v.writelines(my_list)
  v.close()
  print('File for %s has been created' % org) #tells user that function is done

#Obtain list of all annotated experiments, and dictionary where organs are keys for their experiments
header = open('/content/gdrive/My Drive/URECA/ARATH-matrix-filtered.txt','r').readlines()[0] #first row of matrix file
matrix_exps = header.rstrip().split('\t')[1:] #list of exps in matrix
annot = open('/content/gdrive/My Drive/URECA/ARATH-annotation.txt','r').readlines()
exp_list = [] #container for all annotated exps
org_dict = {} #container for organ:exp dictionary

for i in annot:
  cols = i.rstrip().split('\t')
  exp = cols[0] #first word in row is experiment
  org = cols[1].replace(' ', '-') #second word in row is organ, replaces space with hyphen
  if org != '-' and exp in matrix_exps: #ignores exps without annotations and exps not in matrix
    exp_list.append(exp) #adds annotated exps
    org_dict.setdefault(org, set())
    org_dict[org].add(exp) #adds to organ:exp dictionary

#Obtains dictionary where experiments are keys for their tpms
matrix = open('/content/gdrive/My Drive/URECA/ARATH-matrix-filtered.txt','r').readlines()
exp_dict = {} #container for exp:tpm dictionary
ind_list = [] #container for indices of annotated exps
for exp in exp_list:
  ind_list.append(matrix[0].rstrip().split('\t').index(exp))
for line in matrix[1:]: #skips first row, runs through matrix row by row
  for i in range(len(ind_list)): #only takes tpm values of annotated exps
    exp_dict.setdefault(exp_list[i],[])
    exp_dict[exp_list[i]].append(line.rstrip().split('\t')[ind_list[i]]) #adds to exp:tpm dictionary

#Obtains first column of matrix, which is list of genes
gene_col = []
for line in matrix[1:]: #skips first row
  gene = line.split('\t')[0]
  gene_col.append(gene)

#Create file for each organ
for org, exps in org_dict.items():
  file_creator(org, exps)
