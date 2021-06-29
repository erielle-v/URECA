#!/usr/bin/env python3
#This script creates a file with the average TPM values for each gene and organ and a file with SPM values
#Erielle Villanueva - 23 Sep 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/'
import statistics

#Obtain dictionary where organ is key for experiment indices
matrix = open(path + 'ARATH-matrix-filtered.txt','r').readlines()
ind_dict = {} #container for organ:ind dictionary
exps = matrix[0].rstrip().split('\t')[1:] #list of exps in matrix

for i in range(len(exps)):
  org = exps[i].split('+')[-1]
  ind_dict.setdefault(org, set())
  ind_dict[org].add(i+1) #adds to organ:ind dictionary

#Obtain dictionary where organ is key for average tpm values
tpm_dict = {} #container for organ:tpm dictionary
for line in open(path + 'ARATH-matrix-filtered.txt','r'):
  line = line.rstrip().split('\t')
  if line[0] == 'gene': #skips first row
    continue
  for org in ind_dict:
    tpm_dict.setdefault(org, [])
    tpm = statistics.mean([float(line[i]) for i in ind_dict[org]]) #calculates average
    tpm_dict[org].append(tpm)

#Obtain first column of matrix, which is list of genes
gene_col = []
for line in matrix[1:]: #skips first row
  gene = line.split('\t')[0]
  gene_col.append(gene)

#Create file for average TPM values
org_list = sorted(tpm_dict.keys()) #list of organs sorted in alphabetical order
v = open(path + 'Organ-specific-analysis/' + 'ARATH-matrix-average.txt', 'w')
v.write('gene\t' + '\t'.join(org_list) + '\n') #creates first row of the file
for i in range(len(gene_col)):
  gene = gene_col[i]
  tpms = ''
  for org in org_list:
    tpms += '\t' + str(tpm_dict[org][i])
  v.write(gene + tpms + '\n')
v.close()
print('File has been created') #tells user that script is done

#Replace TPM values with SPM values
v = open(path + 'ARATH-matrix-spm.txt', 'w')
for line in open(path + 'ARATH-matrix-average.txt','r'):
  line = line.rstrip().split('\t')
  if line[0] == 'gene':
    v.write('\t'.join(line) + '\n') #creates first row of file
    continue
  gene = line[0]
  total = sum([float(i) for i in line[1:]])
  spm = []
  if total == 0:
    continue #ignores line if the row values are all 0
  else:
    for i in range(1, len(line)):
      spm.append(str(float(line[i]) / total)) #replaces average TPM with SPM
  v.write(gene + '\t' + '\t'.join(spm) + '\n')
v.close()
print('File has been created') #tells user that script is done
