#!/usr/bin/env python

###############################################################################
#                        Plot ROC curve from a CSV file                       #
###############################################################################

# File format required:
# IDs,id1,id2,...
# AUC,auc1,auc2,...
# PRECISION,precision1,precision2,...
# RECALL,recall1,recall2,...
# F-SCORE,f-score1,f-score2,...
# x0,y01,y02,...
# x1,y11,y12,...
# ...

import sys
import matplotlib.pyplot as plt
import csv

if len(sys.argv) != 2:
    sys.exit("Usage: python plot_roc.py roc.csv")

rocFileName = sys.argv[1]
rocFile = open(rocFileName, 'r')
pngFileName = rocFileName + ".png"

print "ROC file: {}".format(rocFileName)
print "PNG file: {}".format(pngFileName)

reader = csv.reader(rocFile, delimiter=",")

ids = next(reader)
aucs = next(reader)
precisions = next(reader)
recalls = next(reader)
fscores = next(reader)
numColumns = len(ids)
numDatasets = numColumns - 1
ids.pop(0)
aucs.pop(0)

columns = []
for i in range(0, numColumns):
    columns.append([])
for row in reader:
    for i in range(0, numColumns):
        columns[i].append(row[i])

plt.figure(1)
for i in range(0, numDatasets):
    plt.plot(columns[0], columns[i+1], label=ids[i]+", AUC="+aucs[i])

plt.title("ROC plot")
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.legend(loc='lower right')
plt.savefig(pngFileName, format='png', dpi=900)

