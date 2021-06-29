#!/usr/bin/env python3
#This script creates a histogram of SPM values (with top 5% threshold) and a file with organ-specific genes
#Erielle Villanueva - 24 Sep 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/Organ-specific-analysis/'
import matplotlib.pyplot as plt

#Obtain a list of all SPM values
spm_list = [] #container for all spm values
for line in open(path + 'ARATH-matrix-spm.txt','r'):
  if line.split('\t')[0] == 'gene': #skips first row
    continue
  spm = line.rstrip().split('\t')[1:]
  spm_list.extend([float(i) for i in spm])

#Output lowest value of top 5%
spm_list.sort(reverse=True) #sorts list in descending order
top_5 = int((0.05*len(spm_list)) - 1)
threshold = spm_list[top_5]
print(threshold)

#Plot histogram of SPM values
plt.hist(spm_list, bins = 50)
plt.xlabel('SPM values', fontsize = 15)
plt.xticks(fontsize = 13)
plt.ylabel('Counts', fontsize = 15)
plt.yticks(fontsize = 13)
plt.axvline(x = threshold, color = 'red')
plt.text(0.41, 10000,'x = %s' % str(threshold), color = 'red', fontsize = 13)
plt.yscale('log')
plt.show()

#Filter genes with spm values >= threshold
v = open(path + 'ARATH-matrix-osg.txt', 'w')
v.write('gene\torgan\n') #creates first row of file
for line in open(path + 'ARATH-matrix-spm.txt','r'):
  gene = line.split('\t')[0]
  if gene == 'gene': #skips first row
    ind_dict = {ind:org for ind, org in enumerate(line.rstrip().split('\t')[1:])} #dictionary of indices to organs
    continue
  spm = line.rstrip().split('\t')[1:]
  osg = []
  for i in range(len(spm)):
    if float(spm[i]) >= threshold: #filters SPM values
      osg.append(ind_dict[i])
  if len(osg) == 0: #adds non-organ-specific genes
    v.write(gene + '\t' + 'Ubiquitous' + '\n')
  elif len(osg) == 1: #adds organ-specific genes
    v.write(gene + '\t' + osg[0] + '\n')
  elif len(osg) > 1: #adds ubiquitous genes
    v.write(gene + '\t' + 'Ubiquitous (' + '/'.join(['%s']*len(osg)) % tuple(osg) + ')\n')
v.close()

print('File has been created') #tells user that script is done
