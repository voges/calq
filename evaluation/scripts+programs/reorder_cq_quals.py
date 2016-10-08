#!/usr/bin/env python

###############################################################################
#           Reorder mapped & unmapped QVs according to a SAM file             #
###############################################################################

import sys

if len(sys.argv) != 4:
    sys.exit("Usage: python reorder_cq_qual.py file.sam mapped.qual unmapped.qual")

samFileName = sys.argv[1]
samFile = open(samFileName, 'r')
mappedQualFileName = sys.argv[2]
mappedQualFile = open(mappedQualFileName, 'r')
unmappedQualFileName = sys.argv[2]
unmappedQualFile = open(unmappedQualFileName, 'r')
newQualFileName = samFileName + '.reordered.qual'
newQualFile = open(newQualFileName, 'w')

print "SAM file: {}".format(samFileName)
print "Mapped QUAL file: {}".format(mappedQualFileName)
print "Unmapped QUAL file: {}".format(unmappedQualFileName)
print "New QUAL file: {}".format(newQualFileName)

idx = 0
headerLines = 0
alignmentLines = 0
mappedLines = 0
unmappedLines = 0

print "Reordering quality values"
while 1:
    print "\r{}".format(idx),
    sys.stdout.flush()

    line = samFile.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        flag = fields[1]
        rname = fields[2]
        pos = fields[3]
        cigar = fields[5]
        seq = fields[9]
        qual = fields[10]

        if (int(flag) & 4) or len(rname) == 0 or rname == "*" or pos == 0 or len(cigar) == 0 or cigar == "*" or len(seq) == 0 or seq == "*" or len(qual) == 0 or qual == "*":
            newQualFile.write(unmappedQualFile.readline())
            unmappedLines += 1
        else:
            newQualFile.write(mappedQualFile.readline())
            mappedLines += 1

        alignmentLines += 1
    else:
        headerLines += 1
    idx += 1
print ""
print "Processed {} lines ({} header lines, {} alignment lines)".format(idx, headerLines, alignmentLines)
print "Of alignment lines, {} were mapped and {} unmapped".format(mappedLines, unmappedLines)

