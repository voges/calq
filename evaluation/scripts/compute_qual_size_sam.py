#!/usr/bin/env python

###############################################################################
#               Compute the quality values size in a SAM file                 #
###############################################################################

import sys

if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.sam".format(sys.argv[0]))

samFileName = sys.argv[1]
samFile = open(samFileName, 'r')

print("SAM file: {}".format(samFileName))

lineCnt = 0
headerLines = 0
alignmentLines = 0
qualityValuesSize = 0

while 1:
    print("\rLine: {}".format(lineCnt)),
    sys.stdout.flush()

    line = samFile.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        qual = fields[10]
        qualityValuesSize += len(qual)
        alignmentLines += 1
    else:
        headerLines += 1
    lineCnt += 1
print("")
print("Processed {} header + {} alignment = {} lines".format(headerLines, alignmentLines, lineCnt))
print("Quality values size: {} bytes".format(qualityValuesSize))
