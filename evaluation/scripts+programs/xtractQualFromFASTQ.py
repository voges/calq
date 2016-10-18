#!/usr/bin/env python

###############################################################################
#               Extract the quality values from a FASTQ file                  #
###############################################################################

import sys


if len(sys.argv) != 2:
    print("Usage: python {} file.fastq 2> file.qual".format(sys.argv[0]))
    sys.exit()

# Get FASTQ file
fastqFileName = sys.argv[1]
if not fastqFileName.endswith(".fastq") and not fastqFileName.endswith(".fq"):
    print("FASTQ file name must end with either '.fastq' or '.fq'")
    sys.exit()
fastqFile = open(fastqFileName, 'r')
print("FASTQ file: {}".format(fastqFileName))

recordCnt = 0
lineCnt = 0

while 1:
    sequenceIdentifierLine = fastqFile.readline().rstrip('\n')
    if not sequenceIdentifierLine:
        break # reached end of file, everything ok
    lineCnt += 1

    rawSequenceLine = fastqFile.readline().rstrip('\n')
    lineCnt += 1

    descriptionLine = fastqFile.readline().rstrip('\n')
    lineCnt += 1

    qualityScoreLine = fastqFile.readline().rstrip('\n')
    lineCnt += 1

    recordCnt += 1

    sys.stderr.write(qualityScoreLine + "\n")

    sys.stdout.flush()

print("Processed {} FASTQ record(s) (i.e. {} lines)".format(recordCnt, lineCnt))

sys.exit()

