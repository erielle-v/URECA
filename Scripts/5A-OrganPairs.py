#!/usr/bin/env python3
#This script creates files for organ-specific genes, connections within organs, and connections between organs
#Erielle Villanueva - 21 May 2021

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/'
from itertools import combinations_with_replacement

#Creates dictionary where organs are keys for their set of specific genes
org_dict = {} #container for org:genes dictionary
for line in open(path + 'Organ-specific-analysis/' + 'ARATH-matrix-osg.txt','r'):
  cols = line.rstrip().split('\t')
  if cols[0] == 'gene': #skips first row
    continue
  if cols[1].startswith('Ubiquitous'):
    continue
  gene, org = cols[0], cols[1]
  org_dict.setdefault(org, set())
  org_dict[org].add(gene)
org_list = sorted(org_dict.keys()) #list of organs sorted in alphabetical order

#Creates files for organ-specific genes (1 per organ, total 9)
for org in org_list:
  u = open(path + 'Functional-analysis/Organ-specific/%s.txt' % org, 'w')
  for gene in org_dict[org]:
    u.write(gene + '\n')
  u.close()

#Creates dictionary for organ pairwise comparisons
org_pairs = sorted(list(combinations_with_replacement(org_list, 2))) #list of all possible organ pairs
pair_dict = {} #container for dictionary where organ pair is key for gene pairs
for pair in org_pairs:
  pair_dict.setdefault(pair, [])

for line in open(path + 'Network-analysis/HRR-NU/' + 'ARATH-matrix-HRRmod-NU.txt', 'r'):
  cols = line.rstrip().split('\t')
  if cols[0] != 'gene_a': #skips first row
    for pair in org_pairs:
      if cols[0] in org_dict[pair[0]] and cols[1] in org_dict[pair[1]]: #if both genes are from the same organ
        pair_dict[pair].append(cols[0] + '\t' + cols[1])
      if pair[0] != pair[1]: #checks other way around if organs are not the same
        if cols[0] in org_dict[pair[1]] and cols[1] in org_dict[pair[0]]:
          pair_dict[pair].append(cols[0] + '\t' + cols[1])

#Creates files for connections within organs (1 per organ, total 9)
for pair in org_pairs:
  if pair[0] == pair[1]: #checks if organs are the same
    v = open(path + 'Functional-analysis/Within-organ/%s.txt' % pair[0], 'w')
    v.write('\n'.join(pair_dict[pair]))
    v.close()

#Creates files for connections between organs (1 per pairwise comparison, total 36)
for pair in org_pairs:
  if pair[0] != pair[1]: #checks if organs are not the same
    w = open(path + 'Functional-analysis/Between-organs/%s-%s.txt' % (pair[0], pair[1]), 'w')
    w.write('\n'.join(pair_dict[pair]))
    w.close()
