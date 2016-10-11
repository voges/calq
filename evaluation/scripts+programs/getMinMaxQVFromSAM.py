#!/usr/bin/env python

###############################################################################
#      Get the minimum and maximum quality values from a SAM file             #
###############################################################################

import sys

if len(sys.argv) != 2:
    sys.exit("Usage: python getMinMaxQVFromSAM.py file.sam")

samFileName = sys.argv[1]
samFile = open(samFileName, 'r')

print "SAM file: {}".format(samFileName)

minQualityScore = sys.maxint
maxQualityScore = -sys.maxint - 1

lineCnt = 0
headerLines = 0
alignmentLines = 0

while 1:
    print "\rLine: {} - qMin:qMax = {}:{}".format(lineCnt, minQualityScore, maxQualityScore),
    sys.stdout.flush()

    line = samFile.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        qual = fields[10]

        if len(qual) > 0:
            for q in qual:
                if ord(q) > maxQualityScore:
                    maxQualityScore = ord(q)
                if ord(q) < minQualityScore:
                    minQualityScore = ord(q)
        else:
            print "Error: No quality scores in line {}".format(lineCnt)

        alignmentLines += 1
    else:
        headerLines += 1
    lineCnt += 1
print ""
print "Processed {} header + {} alignment = {} lines".format(headerLines, alignmentLines, lineCnt)
print "qMin:qMax = {}:{}".format(minQualityScore, maxQualityScore)

