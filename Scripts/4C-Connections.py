#!/usr/bin/env python3
#This script creates a table comparing number of observed and expected connections in organ pairs and a heatmap
#Erielle Villanueva - 23 Mar 2021

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/Network-analysis/'
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

v = open(path + 'ARATH-matrix-organ-pairs-analysis.txt', 'w')
v.write('\t'.join(['Org 1', 'Org 2', 'No. of Observed Connections', 'No. of Expected Connections', 'Ratio of O/E']) + '\n') #header
for line in open(path + 'ARATH-matrix-organ-pairs.txt','r'):
  if line.startswith('Org 1'): #skips first line which is the header
    continue
  line = line.rstrip().split('\t')
  observed = int(line[4]) #number of observed connections
  expected = int(line[2]) * int(line[3]) #number of expected connections
  v.write('\t'.join([line[0], line[1], str(observed), str(expected), str(observed/expected)]) + '\n') #calculates ratio for every row in file
v.close()
print('Organ pairs analysis file has been created') #tells user that organ pairs analysis file is done

#Obtains data for heatmap
data_dict = {} #dict where organ is key for inner dict of organ:ratio
for line in open(path + 'ARATH-matrix-organ-pairs-analysis.txt', 'r'):
  if line.startswith('Org 1'): #skips first line which is the header
    continue
  line = line.rstrip().split('\t')
  org1 = line[0]
  org2 = line[1]
  data_dict.setdefault(org1, {})
  data_dict.setdefault(org2, {})
  data_dict[org1].update({org2:line[4]})
  data_dict[org2].update({org1:line[4]})
org_list = sorted(data_dict.keys())

heatmap_data = [] #container for 2D list of values
for org1 in org_list:
  temp_list = []
  for org2 in data_dict[org1]:
    temp_list.append(float(data_dict[org1][org2]) * 1000) #normalise values
  heatmap_data.append(temp_list)

#Edits list of organs for axes_labels
for i in range(len(org_list)):
  if org_list[i].endswith('meristem'):
    org_list[i] = org_list[i].split('-')[0] + ' meristem' #replaces - with space

#Creates heatmap
mask = np.ones_like(heatmap_data, dtype=np.bool)
mask[np.tril_indices_from(mask)] = False #only shows lower part below diagonal
color = sns.color_palette("Blues", as_cmap = True)
sns.heatmap(heatmap_data, vmax = 0.2, xticklabels = org_list, yticklabels = org_list, cmap = color, annot = True, mask = mask)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.text(9.5, 9.4, '$\mathregular{x10^{-3}}$')
plt.rcParams.update({'font.size': 12})
plt.show()
