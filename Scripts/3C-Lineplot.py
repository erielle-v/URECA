#!/usr/bin/env python3
#This script creates a line plot of percentages of highly/moderately co-expressed genes for HRR/MR
#Erielle Villanueva - Oct 15 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/Network-analysis/'
import matplotlib.pyplot as plt

#Creates list of high and mod genes for HRR
hrr_high = []
hrr_mod = []
org_list = []
for line in open(path + 'HRR/' + 'ARATH-matrix-HRRper.txt','r'):
  cols = line.rstrip().split('\t')
  if cols[0] == 'organ' or cols[0] == 'Ubiquitous':
    continue
  org, high, mod = cols[0], cols[3], cols[4]
  org_list.append(org.replace('-', ' '))
  hrr_high.append(float(high))
  hrr_mod.append(float(mod))

#Creates list of high and mod genes for MR
mr_high = []
mr_mod = []
for line in open(path + 'MR/' + 'ARATH-matrix-MRper.txt','r'):
  cols = line.rstrip().split('\t')
  if cols[0] == 'organ' or cols[0] == 'Ubiquitous':
    continue
  high, mod = cols[3], cols[4]
  mr_high.append(float(high))
  mr_mod.append(float(mod))

#Creates scatter plot
plt.plot(org_list, hrr_high, label = 'HRR High')
plt.plot(org_list, hrr_mod, label = 'HRR Mod')
plt.plot(org_list, mr_high, label = 'MR High')
plt.plot(org_list, mr_mod, label = 'MR Mod')
plt.xticks(rotation = 50, ha = "right", fontsize = 15)
plt.yticks(fontsize = 13)
plt.ylabel('% of genes out of total', fontsize = 15)
plt.legend(loc = 'center', ncol = 2, bbox_to_anchor=(0.48, 1.15), fontsize = 15)
fig = plt.gcf()
fig.set_size_inches(7, 4)
plt.show()
