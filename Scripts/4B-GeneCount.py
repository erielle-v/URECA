#!/usr/bin/env python3
#This script creates a table comparing number of connected genes in organ pairs and number of genes in each organ
#Erielle Villanueva - 14 Mar 2021

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/'
from itertools import combinations_with_replacement

#This function takes the input of the names of 2 organs as strings, and outputs a string with the number of genes in each organ
def gene_counter(org1, org2):
  count = 0 #counter for number of connected genes 
  for pair in pair_set: #checks every pair to see if they belong to (org1, org2) or (org2, org1)
    if pair[0] in org_dict[org1] and pair[1] in org_dict[org2]:
      count += 1
      if org1 == org2: #if the organs are the same, skips to the next gene pair to avoid double counting
        continue
    if pair[0] in org_dict[org2] and pair[1] in org_dict[org1]:
      count += 1
  org1_len = str(len(org_dict[org1])) #number of genes in org1
  org2_len = str(len(org_dict[org2])) #number of genes in org1
  row = '\t'.join( [org1, org2, org1_len, org2_len, str(count)] ) + '\n' #row that will be written in the file
  return row

#Obtains dictionary for organs and their organ-specific genes
org_dict = {} #container for org:genes dictionary
for line in open(path + 'Organ-specific-analysis/' + 'ARATH-matrix-osg.txt','r'):
  cols = line.rstrip().split('\t')
  gene, org = cols[0], cols[1]
  if gene != 'gene' and org.startswith('Ubiquitous') == False: #skips first row and ignores ubiquitous genes
    org_dict.setdefault(org, set())
    org_dict[org].add(gene)
org_list = sorted(org_dict.keys()) #list of organs sorted in order, to order the organs in the files consistently

pair_set = set() #container for set of all gene pairs in network
for line in open(path + 'Network-analysis/HRR-NU/' + 'ARATH-matrix-HRRmod-NU.txt', 'r'):
  cols = line.rstrip().split('\t')
  if cols[0] != 'gene_a': #skips first row
    gene_pair = (cols[0], cols[1])
    pair_set.add(gene_pair)

v = open(path + 'Network-analysis/' + 'ARATH-matrix-organ-pairs.txt', 'w')
v.write('\t'.join(['Org 1', 'Org 2', 'No. of Org 1 Genes', 'No. of Org 2 Genes', 'No. of Connections']) + '\n') #header
org_pairs = list(combinations_with_replacement(org_list, 2)) #list of all possible organ pairs
for pair in org_pairs:
  v.write(gene_counter(pair[0], pair[1]))
v.close()
print('Organ pairs file has been created') #tells user that organ pairs file is done
