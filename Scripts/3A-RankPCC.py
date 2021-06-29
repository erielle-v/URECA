#!/usr/bin/env python3
#This script calculates the HRR/MR (PCC ranking methods) of gene pair PCCs
#Erielle Villanueva - 29 Sep 2020

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

#Creates set of organ-specific genes
osg_set = set()
for line in open(path + 'Organ-specific-analysis/' + 'ARATH-matrix-osg.txt','r'):
  cols = line.rstrip().split('\t')
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
v = open(path + 'Network-analysis/HRR/' + 'ARATH-matrix-HRR.txt','w')
for line in open(path + 'Organ-PCC/' + 'ARATH-matrix-PCC.txt','r'):
  gene1, gene2, r = line.rstrip().split('\t')
  if gene1 in osg_set and gene2 in osg_set: #only adds if genes are in organ-specific file
    rank = str(max(master_dict[gene1][gene2],master_dict[gene2][gene1])) #calculates HRR
    v.write('\t'.join([gene1, gene2, rank, r]) + '\n')
v.close()
print('HRR File has been created') #tells user that HRR file is done

#Create MR file
v = open(path + 'Network-analysis/MR/' + 'ARATH-matrix-MR.txt','w')
for line in open(path + 'Organ-PCC/' + 'ARATH-matrix-PCC.txt','r'):
  gene1, gene2, r = line.rstrip().split('\t')
  if gene1 in osg_set and gene2 in osg_set: #only adds if genes are in organ-specific file
    rank = str(math.sqrt(master_dict[gene1][gene2]*master_dict[gene2][gene1])) #calculates MR
    v.write('\t'.join([gene1, gene2, rank, r]) + '\n')
v.close()
print('MR File has been created') #tells user that MR file is done
