#!/usr/bin/env python3
#This script combines Prof Marek's script and my script to output filtered PCCs for organ files and complete matrix file
#Erielle Villanueva - 15 Sep 2020

#Mount Google Drive
from google.colab import drive
drive.mount('/content/gdrive')
import numpy as np
import sys
patha = '/content/gdrive/My Drive/URECA/'

def pcc_file_creator(filename):
    output = path + 'Organ-PCC/' + '%s-PCC.txt' % filename.split('.')[0]

    # Read Matrix and store nominators and denominators
    with open(path + '/Organ' + '%s.txt' % filename, 'r') as fin:
        nominators, denominators, genes = [], [], []

        header = fin.readline()
        size = len(header.strip().split('\t'))

        for line in fin:
                parts = line.rstrip().split("\t")

                if size != len(parts):
                    print("Warning! Unequal number of columns found in line:\n%s.\nExpression matrix corrupt. Aborting!\n" % line, file=sys.stderr)
                    quit()

                if len(parts) == size:
                    temp = []
                    for j in range(1, len(parts)):
                        try:
                            temp.append(float(parts[j]))
                        except ValueError:
                            print("Warning! Non-number character found in line:\n%s.\nExpression matrix corrupt. Aborting!\n" % line, file=sys.stderr)
                            quit()

                    row_values = np.array(temp)
                    nomi = row_values-(sum(row_values)/len(row_values))
                    denomi = np.sqrt(sum(nomi**2))

                    if denomi != 0.0:
                        nominators.append(nomi)
                        denominators.append(denomi)
                        genes.append(parts[0])

    nominators = np.array(nominators)
    denominators = np.array(denominators)

    # Calculate PCC and write output
    with open(output, 'w') as fout: #, open(mcl_output, 'w') as mcl_out
        print("Database OK.\nCalculating Pearson Correlation Coefficient and ranks.\n")
        for i, (nom, denom, gene) in enumerate(zip(nominators, denominators, genes), start=1):

            nominator = np.dot(nominators, nom)
            denominator = np.dot(denominators, denom)
            pcc_values = nominator/denominator

            data = [{'score': p,
                    'gene': g,
                    'string': g + '\t' + str(p)} for g, p in zip(genes, pcc_values) if g != gene and p > 0.8]

            for s in data:
                if gene > s['string']:
                    fout.writelines(gene + '\t' + s['string'] + "\n")

    print("PCCs calculated and saved as %s." % (output), file=sys.stderr)

#Create PCC file for each organ and for complete matrix
filelist = ['Apical-meristem', 'Female', 'Flower', 'Leaf', 'Male', \
            'Root-meristem', 'Root', 'Seeds', 'Stem', 'ARATH-matrix-filtered']
for filename in filelist:
  pcc_file_creator(filename)
