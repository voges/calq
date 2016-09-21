import sys


if len(sys.argv) != 2:
    sys.exit("Usage: python validate_fastq.py file.fastq")

# Get FASTQ file
fastqFileName = sys.argv[1]
if not fastqFileName.endswith(".fastq") and not fastqFileName.endswith(".fq"):
    sys.exit("FASTQ file name must end with either '.fastq' or '.fq'")
fastqFile = open(fastqFileName, 'r')
print "FASTQ file: {}".format(fastqFileName)

minQualityScore = sys.maxint
maxQualityScore = -sys.maxint - 1

recordCnt = 0
lineCnt = 0

while 1:
    # Try to read a whole record (i.e. 4 lines)
    sequenceIdentifierLine = fastqFile.readline().rstrip('\n')
    if not sequenceIdentifierLine:
        break  # reached end of file, everything ok
    lineCnt += 1

    rawSequenceLine = fastqFile.readline().rstrip('\n')
    if not rawSequenceLine:
        print "Parse error (record {}, line {}): Record is not complete".format(recordCnt+1, lineCnt+1)
        break
    lineCnt += 1

    descriptionLine = fastqFile.readline().rstrip('\n')
    if not descriptionLine:
        print "Parse error (record {}, line {}): Record is not complete".format(recordCnt + 1, lineCnt + 1)
        break
    lineCnt += 1

    qualityScoreLine = fastqFile.readline().rstrip('\n')
    if not qualityScoreLine:
        print "Parse error (record {}, line {}): Record is not complete".format(recordCnt + 1, lineCnt + 1)
        break
    lineCnt += 1

    recordCnt += 1

    # Check quality score line

    # Typical ranges:
    #  Sanger         Phred+33   [0,40]
    #  Solexa         Solexa+64  [-5,40]
    #  Illumina 1.3+  Phred+64   [0,40]
    #  Illumina 1.5+  Phred+64   [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
    #  Illumina 1.8+  Phred+33   [0,41]

    if len(qualityScoreLine) > 0:
        for q in qualityScoreLine:
            if ord(q) > maxQualityScore:
                maxQualityScore = ord(q)
            if ord(q) < minQualityScore:
                minQualityScore = ord(q)
    else:
        print "Error (record {}, line {}): Quality score line is empty".format(recordCnt, lineCnt)

    print "\rProcessed {} record(s) - qMin:qMax = {}:{}".format(recordCnt, minQualityScore, maxQualityScore),
    sys.stdout.flush()

print "Processed {} FASTQ record(s) (i.e. {} lines)".format(recordCnt, lineCnt)
print "qMin:qMax = {}:{}".format(minQualityScore, maxQualityScore)

sys.exit()
