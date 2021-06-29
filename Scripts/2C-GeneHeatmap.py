#!/usr/bin/env python3
#This script creates a file with normalised TPM values and a heat map of TPM values for each organ
#Erielle Villanueva - 22 Oct 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/Organ-specific-analysis/'
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#Creates dictionary where organs are keys for their set of specific genes
org_dict = {} #container for org:genes dictionary
for line in open(path + 'ARATH-matrix-osg.txt','r'):
  cols = line.rstrip().split('\t')
  if cols[0] == 'gene': #skips first row
    continue
  if cols[1].startswith('Ubiquitous'): #filters for genes with >1 organ
    continue
  gene, org = cols[0], cols[1]
  org_dict.setdefault(org, set())
  org_dict[org].add(gene)
org_list = sorted(org_dict.keys()) #list of organs sorted in alphabetical order

#Creates file with average TPM values for organ-specific genes
v = open(path + 'osg-average.txt', 'w')
for org in org_list:
  for line in open(path + 'ARATH-matrix-average.txt','r'):
    gene = line.rstrip().split('\t')[0]
    if gene == 'gene': #skips first row
      v.write(line) #creates first row of file
      continue
    if gene in org_dict[org]:
      v.write(line)
v.close()
print('Osg average TPM file has been created') #tells user that osg average TPM file is done

#Replaces TPM values with normalised TPM values
v = open(path + 'ARATH-matrix-norm.txt', 'w')
for line in open(path + 'osg-average.txt','r'):
  line = line.rstrip().split('\t')
  gene = line[0]
  if gene == 'gene':
    v.write('\t'.join(line) + '\n') #creates first row of file
    continue
  highest = max([float(i) for i in line[1:]])
  norm = []
  if highest == 0:
    continue #ignores line if the row values are all 0
  else:
    for i in range(1, len(line)):
      norm.append(str(float(line[i]) / highest)) #replaces average TPM with normalised TPM
  v.write(gene + '\t' + '\t'.join(norm) + '\n')
v.close()
print('Normalised TPM file has been created') #tells user that normalised TPM file is done

#Obtains data for heatmap
heatmap_data = [] #container for 2D list of normalised TPM values
for line in open(path + 'ARATH-matrix-norm.txt', 'r'):
  line = line.rstrip().split('\t')
  if line[0] == 'gene':
    org_list = line[1:]
    continue
  data = [float(i) for i in line[1:]] #obtains the TPM values for each gene
  heatmap_data.append(data)

#Edits list of organs for x-axis labels
for i in range(len(org_list)):
  count = len(org_dict[org_list[i]])
  if org_list[i].endswith('meristem'):
    org_list[i] = org_list[i].split('-')[0] + ' meristem'
  org_list[i] += ' (%s)' % count

#Creates heatmap
color = sns.color_palette("Blues", as_cmap = True)
heatmap = sns.heatmap(heatmap_data, xticklabels = org_list, cmap = color, cbar_kws = dict(use_gridspec=False,location="top"))
plt.xticks(rotation = 50, ha = "right")
plt.yticks(np.arange(len(heatmap_data), step = 1000), [i for i in range(0, len(heatmap_data), 1000)])
plt.rcParams.update({'font.size': 15})
fig = plt.gcf()
fig.set_size_inches(7, 7)
plt.show()
