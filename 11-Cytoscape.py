#!/usr/bin/env python3
#This script creates HRR files with mod and high gene pairs and a features files with mod genes and their organs
#Final version, removes 'ubiquitous' genes and uses HRR cut-offs of 2 and 9 for mod and high respectively
#Erielle Villanueva - 04 Nov 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/'
import math

#Function that outputs gene:rank dictionary, takes gene:pcc dictionary as argument
def rank(gene_dict):
  pcc = sorted(gene_dict, key=gene_dict.get, reverse=True) #creates list of genes sorted by PCCs, highest to lowest
  rank_dict = {gene:rank for rank, gene in enumerate(pcc, 1)} #dictionary where gene is key for their rank
  return rank_dict

#Creates set of organ-specific genes without 'ubiquitous' genes
osg_set = set()
for line in open(path + 'Organ-specific-analysis/' + 'ARATH-matrix-osg.txt','r'):
  cols = line.rstrip().split('\t')
  if cols[1].startswith('Ubiquitous'): #filters for 'ubiquitous' genes
    continue
  gene = cols[0]
  osg_set.add(gene)

#Obtain dictionary where gene1 is key for dictionary of gene2:pcc, and vice versa
master_dict = {}
for line in open(path + 'Organ-PCC/' + 'ARATH-matrix-PCC.txt','r'):
  gene1, gene2, r = line.rstrip().split('\t')
  if gene1 in osg_set and gene2 in osg_set: #only adds if genes are in organ-specific file
    master_dict.setdefault(gene1,{})
    master_dict.setdefault(gene2,{})
    master_dict[gene1].update({gene2:r})
    master_dict[gene2].update({gene1:r})
for i in master_dict: #replaces gene:pcc dictionary with gene:rank dictionary
  master_dict[i] = rank(master_dict[i])

#Create HRR file
v = open(path + 'Network-analysis/HRR-NU/' + 'ARATH-matrix-HRR-NU.txt','w')
for line in open(path + 'Organ-PCC/' + 'ARATH-matrix-PCC.txt','r'):
  gene1, gene2, r = line.rstrip().split('\t')
  if gene1 in osg_set and gene2 in osg_set: #only adds if genes are in organ-specific file
    rank = str(max(master_dict[gene1][gene2],master_dict[gene2][gene1])) #calculates HRR
    v.write('\t'.join([gene1, gene2, rank, r]) + '\n')
v.close()
print('HRR-NU File has been created') #tells user that HRR file is done

#Creates dictionary where organs are keys for their set of specific genes
org_dict = {} #container for org:genes dictionary
for line in open(path + 'Organ-specific-analysis/' + 'ARATH-matrix-osg.txt','r'):
  cols = line.rstrip().split('\t')
  if cols[0] == 'gene': #skips first row
    continue
  if cols[1].startswith('Ubiquitous'): #filters for 'ubiquitous' genes
    continue
  gene, org = cols[0], cols[1]
  org_dict.setdefault(org, set())
  org_dict[org].add(gene)

#Sorts values in HRR file and finds lowest value of top 1% and 5%
hrr_list = [] #container for list of HRR values
for line in open(path + 'Network-analysis/HRR/' + 'ARATH-matrix-HRR.txt','r'):
  hrr = line.rstrip().split('\t')[2]
  hrr_list.append(int(hrr))
hrr_list.sort() #sorts list in ascending order
top_5 = int((0.05*len(hrr_list)) - 1) #lowest value of 5%
top_1 = int((0.01*len(hrr_list)) - 1) #lowest value of 1%

#Creates HRR files
u = open(path + 'Network-analysis/HRR-NU/' + 'ARATH-matrix-HRRhigh-NU.txt','w')
u.write('\t'.join(['gene_a', 'gene_b', 'hrr', 'pcc']) + '\n')
v = open(path + 'Network-analysis/HRR-NU/' + 'ARATH-matrix-HRRmod-NU.txt','w') #HRR mod file
v.write('\t'.join(['gene_a', 'gene_b', 'hrr', 'pcc']) + '\n')
w = open(path + 'Network-analysis/HRR-NU/' + 'ARATH-matrix-HRRfeatures-NU.txt','w') #HRR features file
w.write('gene\torgan\n')
hrr_set = set() #container for set of HRR mod genes
for line in open(path + 'Network-analysis/HRR-NU/' + 'ARATH-matrix-HRR-NU.txt','r'):
  g1, g2, hrr, r = line.rstrip().split('\t')
  if hrr_list[top_1] >= int(hrr): #filters highly co-expressed genes
    u.write(line)
  if hrr_list[top_5] >= int(hrr): #filters moderately co-expressed genes
    v.write(line)
    hrr_set.add(g1)
    hrr_set.add(g2)
u.close()
v.close()
print('HRR-NU mod and high files have been created') #tells user that HRR mod and high files are done
for gene in hrr_set:
  for org in org_dict:
    if gene in org_dict[org]:
      w.write(gene + '\t' + org + '\n')
w.close()
print('HRR-NU features file has been created') #tells user that HRR features file is done
