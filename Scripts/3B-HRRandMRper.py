#!/usr/bin/env python3
#This script checks what percentage of organ-specific genes are highly/moderately co-expressed for HRR and MR
#Erielle Villanueva - Oct 13 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/'

#Creates dictionary where organs are keys for their set of genes
org_dict = {} #container for org:genes dictionary
for line in open(path + 'Organ-specific-analysis/' + 'ARATH-matrix-osg.txt','r'):
  cols = line.rstrip().split('\t')
  if cols[0] == 'gene': #skips first row
    continue
  if cols[1].startswith('Ubiquitous') and len(cols[1]) > 10: #filters for genes with >1 organ
    cols[1] = 'Ubiquitous'
  gene, org = cols[0], cols[1]
  org_dict.setdefault(org, set())
  org_dict[org].add(gene)
org_list = sorted(org_dict.keys()) #list of organs sorted in alphabetical order

#Obtains dictionaries for HRR moderate and high genes
mod_dict = {org:0 for org in org_list} #container for dictionary where organs are keys for number of mod genes
mod_set = set() #container for set of mod geens
for line in open(path + 'Network-analysis/HRR/' + 'ARATH-matrix-HRRmod.txt'):
  if line.startswith('gene_a'):
    continue
  g1, g2, hrr, r = line.rstrip().split('\t')
  mod_set.add(g1)
  mod_set.add(g2)
for gene in mod_set:
  for org in org_dict:
    if gene in org_dict[org]:
      mod_dict[org] += 1

high_dict = {org:0 for org in org_list} #container where organs are keys for number of high genes
high_set = set() #container for set of mod genes
for line in open(path + 'Network-analysis/HRR/' + 'ARATH-matrix-HRRhigh.txt'):
  if line.startswith('gene_a'):
    continue
  g1, g2, hrr, r = line.rstrip().split('\t')
  high_set.add(g1)
  high_set.add(g2)
for gene in high_set:
  for org in org_dict:
    if gene in org_dict[org]:
      high_dict[org] += 1

#Creates HRR percentage file
v = open(path + 'Network-analysis/HRR/' + 'ARATH-matrix-HRRper.txt','w')
v.write('\t'.join(['organ','no. of genes high', 'no. of genes mod', '%high', '%mod']) + '\n')
for org in org_list:
  num_high = high_dict[org]
  num_mod = mod_dict[org]
  per_high = (num_high / len(org_dict[org])) * 100
  per_mod = (num_mod / len(org_dict[org])) * 100
  v.write('\t'.join([org, str(num_high), str(num_mod), str(per_high), str(per_mod)]) + '\n')
v.close()
print('HRR percentage file has been created') #tells user that HRR percentage file is done

#Obtains dictionaries for MR moderate and high genes
mod_dict = {org:0 for org in org_list} #container for dictionary where organs are keys for number of mod genes
mod_set = set() #container for set of mod geens
for line in open(path + 'Network-analysis/MR/' + 'ARATH-matrix-MRmod.txt'):
  if line.startswith('gene_a'):
    continue
  g1, g2, mr, r = line.rstrip().split('\t')
  mod_set.add(g1)
  mod_set.add(g2)
for gene in mod_set:
  for org in org_dict:
    if gene in org_dict[org]:
      mod_dict[org] += 1

high_dict = {org:0 for org in org_list} #container where organs are keys for number of high genes
high_set = set() #container for set of mod geens
for line in open(path + 'Network-analysis/MR/' + 'ARATH-matrix-MRhigh.txt'):
  if line.startswith('gene_a'):
    continue
  g1, g2, mr, r = line.rstrip().split('\t')
  high_set.add(g1)
  high_set.add(g2)
for gene in high_set:
  for org in org_dict:
    if gene in org_dict[org]:
      high_dict[org] += 1

#Creates MR percentage file
v = open(path + 'Network-analysis/MR/' + 'ARATH-matrix-MRper.txt','w')
v.write('\t'.join(['organ','no. of genes high', 'no. of genes mod', '%high', '%mod']) + '\n')
for org in org_list:
  num_high = high_dict[org]
  num_mod = mod_dict[org]
  per_high = (num_high / len(org_dict[org])) * 100
  per_mod = (num_mod / len(org_dict[org])) * 100
  v.write('\t'.join([org, str(num_high), str(num_mod), str(per_high), str(per_mod)]) + '\n')
v.close()
print('MR percentage file has been created') #tells user that MR percentage file is done
