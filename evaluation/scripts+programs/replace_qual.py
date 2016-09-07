#!/usr/bin/env python

###############################################################################
#                    Replace quality values in a SAM file                     #
###############################################################################

import sys

if len(sys.argv) != 3:
    sys.exit("Usage: python replace_qual.py file.sam new_qual.qual")

samFileName = sys.argv[1]
samFile = open(samFileName, 'r')
newQualFileName = sys.argv[2]
newQualFile = open(newQualFileName, 'r')
newSAMFileName = samFileName + '.new_qual.sam'
newSAMFile = open(newSAMFileName, 'w')

print "SAM file: {}".format(samFileName)
print "New QUAL file: {}".format(newQualFileName)
print "New SAM file: {}".format(newSAMFileName)

idx = 0
headerLines = 0
alignmentLines = 0

print "Replacing quality values"
while 1:
    print "\r{}".format(idx),
    sys.stdout.flush()

    line = samFile.readline()
    if not line:
        break
    if line[0] != '@':
        newQual = newQualFile.readline()
        a = line.split('\t')
        newLine = a[0] + '\t' + a[1] + '\t' + a[2] + '\t' + a[3] + '\t' + a[4] + '\t' + a[5] + '\t'
        newLine += a[6] + '\t' + a[7] + '\t' + a[8] + '\t' + a[9] + '\t'
        newQual = newQual.rstrip('\n')
        newLine += newQual
        if len(a) > 11:
            b = 12
            l = a[11]
            while b < len(a):
                l += '\t' + a[b]
                b += 1
            newLine += '\t' + l
        else:
            newLine += '\n'
        newSAMFile.write(newLine)
        alignmentLines += 1
    else:
        newSAMFile.write(line)
        headerLines += 1
    idx += 1
print ""
print "Processed {} lines ({} header lines, {} alignment lines)".format(idx, headerLines, alignmentLines)
