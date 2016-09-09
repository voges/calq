#!/usr/bin/env python

###############################################################################
#                            Plot calq statistics                             #
###############################################################################

# File format required:
# rname,pos,depth,entropy,k
# rn1,0,3,0.1,3

import sys
import matplotlib.pyplot as plt
import csv

if len(sys.argv) != 2:
    sys.exit("Usage: python plot_calq_stats.py stats.cq.stats")

statsFileName = sys.argv[1]
statsFile = open(statsFileName, 'r')
pngFileName = statsFileName + ".png"
epsFileName = statsFileName + ".eps"

print "STATS file: {}".format(statsFile)
print "PNG file: {}".format(pngFileName)
print "EPS file: {}".format(epsFileName)

reader = csv.reader(statsFile, delimiter=",")

ids = next(reader)
aucs = next(reader)
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
plots = []
for i in range(0, numDatasets):
    plt.plot(columns[0], columns[i+1], label=ids[i]+", AUC="+aucs[i])

plt.title("ROC plot")
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.legend(loc='lower right')
plt.savefig(pngFileName, format='png')
plt.savefig(epsFileName, format='eps')
