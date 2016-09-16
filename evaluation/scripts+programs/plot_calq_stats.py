#!/usr/bin/env python

###############################################################################
#                       Plot calq encoder statistics                          #
###############################################################################

# File format required:
# rname,pos,depth,entropy,k
# rn1,0,3,0.1,3
# ...

import sys
import matplotlib.pyplot as plt
import csv
from math import log10

if len(sys.argv) != 2:
    sys.exit("Usage: python plot_calq_stats.py stats.cq.stats")

statsFileName = sys.argv[1]
statsFile = open(statsFileName, 'r')
pngFileName = statsFileName + ".png"

print "STATS file: {}".format(statsFile)
print "PNG file: {}".format(pngFileName)

reader = csv.reader(statsFile, delimiter=",")

headers = next(reader)
numColumns = len(headers)

columns = []
for i in range(0, numColumns):
    columns.append([])

for row in reader:
    for i in range(0, numColumns):
        columns[i].append(row[i])

plt.figure(1)

plt.plot(columns[1], columns[2], label=headers[2])

tmp = []
for elem in columns[3]:
    tmp.append(-log10(float(elem)))

plt.plot(columns[1], tmp, label=headers[3])

plt.title("Calq encoder statistics")
plt.xlabel('Locus')
plt.ylabel('Entropy')
plt.legend(loc='lower right')
plt.savefig(pngFileName, format='png', dpi=900)

