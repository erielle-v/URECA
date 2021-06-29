#!/usr/bin/env python3
#This script creates a clustered heat map and a filtered complete matrix file
#Erielle Villanueva - 03 Dec 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
path = '/content/gdrive/My Drive/URECA/'
import pandas as pd
from sklearn.preprocessing import StandardScaler
import seaborn as sns

#Converts matrix to dataframe
df = pd.read_csv(path + 'ARATH-matrix-labelled.txt', sep = '\t', header = 0)
tpm = df.iloc[:,1:] #selects only the tpm values
scaler = StandardScaler()
scaler.fit(tpm) 
df2 = pd.DataFrame(scaler.transform(tpm), columns = tpm.columns.values) #standardises the data
df3 = df2.corr() #calculates the pcc between every column

#Filters matrix by PCC
nodrop_set = set() #container for set of samples that pass the filter
for col in df3: #iterates through every column
  series = df3[col].sort_values(ascending = False) #sorts each column according to pcc
  org1 = series.name.split(', ')[1]
  for row in series.index: #iterates through every row
    pcc = series[row] 
    org2 = row.split(', ')[1]
    if col != row: #does not compare sample to itself
      if pcc >= 0.95 and org2 == org1: #filters samples if PCC >= 0.95 and same organ
        nodrop_set.add(row)
drop_list = list(set(series.index) - nodrop_set) #obtains list of samples that did not pass the filter
df2.drop(columns = drop_list, inplace = True) #removes samples from dataframe
df4 = df2.corr() #re-calculates pcc between every column in new dataframe

#Creates list of colors for color bar
col_dict = {'Apical meristem':'#ff3333', 'Female':'#ff9900', 'Flower':'#ffff66', 'Leaf':'#66ff66', 'Male': '#33ffff',\
            'Root':'#0000ff', 'Root meristem':'#9900ff', 'Seeds':'#ff66ff', 'Stem':'#993300'} #dictionary of org:colour
col_list = [] #container for list of colours
for sample in df4.columns:
  org = sample.split(', ')[1]
  for j in col_dict:
    if org == j:
      col_list.append(col_dict[j]) #adds colours to list according to index

#Plots clustered heat map
color = sns.color_palette("ch: s = .25, rot = -.2", as_cmap = True) #sets colors of heatmap
heatmap = sns.clustermap(df4, cmap = color, cbar_pos = (0.02, 0.8, 0.03, 0.18), xticklabels = False, figsize = (8.5, 6.5), row_colors = col_list) #plots map
heatmap.ax_col_dendrogram.set_visible(False) #removes x-axis dendrogram
heatmap.ax_heatmap.set_title('Arabidopsis thaliana', fontsize = 20, fontstyle = 'italic') #adds title
for label in col_dict: 
  heatmap.ax_row_dendrogram.bar(0, 0, color=col_dict[label], label=label, linewidth=0)
heatmap.ax_row_dendrogram.legend(loc = 'center', ncol = 3, bbox_to_anchor=(3.15, 1.18), fontsize = 12) #sets legend for color bar

#Obtains list of indices for samples to be removed
header = open(path + 'ARATH-matrix-labelled.txt', 'r').readlines()[0] #reads first line of file
samples = header.rstrip().split('\t')
ind_list = [] #container for indices of samples that passed the check
for i in range(len(samples)):
  if samples[i] not in drop_list:
    ind_list.append(i)

#Creates new filtered matrix file
v = open(path + 'ARATH-matrix-filtered.txt', 'w')
for line in open(path + 'ARATH-matrix-labelled.txt', 'r'): #goes through every row in labelled matrix file
  tpms = line.rstrip().split('\t')
  write_list = []
  for i in ind_list: #only adds columns of samples that passed the check
    write_list.append(tpms[i])
  v.write('\t'.join(write_list) + '\n')
v.close()
print('Filtered matrix file has been created') #tells user that matrix file is done
