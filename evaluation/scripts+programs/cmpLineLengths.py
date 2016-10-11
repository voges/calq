#!/usr/bin/env python

###############################################################################
#                    Compare the line lengths in two files                    #
###############################################################################

import sys

if len(sys.argv) != 3:
    sys.exit("Usage: python cmpLineLengths.py file1 file2")

firstFileName = sys.argv[1]
firstFile = open(firstFileName, 'r')
secondFileName = sys.argv[2]
secondFile = open(secondFileName, 'r')

print "1st file: {}".format(firstFileName)
print "2nd file: {}".format(secondFileName)

lineCnt = 0

print "Comparing line lengths"
while 1:
    print "\r{}".format(lineCnt),
    sys.stdout.flush()

    line1 = firstFile.readline()
    line2 = secondFile.readline()
    line1 = line1.rstrip('\n')
    line2 = line2.rstrip('\n')

    if len(line1) != len(line2):
        print ""
        print "Lines {} do not have same length".format(lineCnt)
        print "1 (len={}): {}".format(len(line1), line1)
        print "2 (len={}): {}".format(len(line2), line2)
        break

    lineCnt += 1
print ""
print "Processed {} lines".format(lineCnt)

