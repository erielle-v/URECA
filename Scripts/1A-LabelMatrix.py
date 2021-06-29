#!/usr/bin/env python3
#This script adds organ annotations to sample IDs in the complete matrix
#Erielle Villanueva - 26 Nov 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/'

#Create dictionary where exp IDs are keys for their organs
exp_dict = {} #container for exp:organ dictionary
missing = ['SRR3993761', 'SRR1005241', 'SRR1005238', 'SRR1005240', 'SRR1005239']
for line in open(path + 'ARATH-annotation.txt', 'r'):
  cols = line.rstrip().split('\t')
  exp = cols[0] #first word in row is experiment
  org = cols[1].replace(' ', '-') #second word in row is organ, replaces space with hyphen
  if org != '-' and exp not in missing: #ignores exps without annotations and missing exps
    exp_dict[exp] = org #adds to exp:organ dictionary

#Creates copy of matrix file
v = open(path + 'ARATH-matrix-labelled.txt', 'w')
ind_list = [] #container for list of indices of columns with annotated organs
for line in open(path + 'ARATH-matrix.txt', 'r'):
  exp = line.rstrip().split('\t')
  if exp[0] == 'gene':
    headers = ['gene'] #first row of file
    for i in range(len(exp[1:])):
      if exp[i+1] in exp_dict:
        headers.append(exp[i+1] + ', ' + exp_dict[exp[i+1]].replace('-', ' '))
        ind_list.append(i+1) 
    v.write('\t'.join(headers) + '\n')
    continue
  row = exp[0]
  for ind in ind_list: #writes columns with annotated organs
    row += '\t' + exp[ind]
  row += '\n'
  v.write(row)
v.close()
print('Copy of matrix file has been created') #tells user that matrix file is done
