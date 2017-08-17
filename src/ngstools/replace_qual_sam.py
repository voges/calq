#!/usr/bin/env python

import sys

# Usage
if len(sys.argv) != 3:
    sys.exit("Usage: python {} file.sam new.qual".format(sys.argv[0]))

sam_file_name = sys.argv[1]
sam_file = open(sam_file_name, 'r')
new_qual_file_name = sys.argv[2]
new_qual_file = open(new_qual_file_name, 'r')

while 1:
    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        new_qual = new_qual_file.readline()
        a = line.split('\t')
        new_line = a[0] + '\t' + a[1] + '\t' + a[2] + '\t' + a[3] + '\t' + a[4] + '\t' + a[5] + '\t'
        new_line += a[6] + '\t' + a[7] + '\t' + a[8] + '\t' + a[9] + '\t'
        new_qual = new_qual.rstrip('\n')
        new_line += new_qual
        if len(a) > 11:
            for i in range(11, len(a)):
                new_line += '\t' + a[i]
        else:
            new_line += '\n'
        sys.stdout.write(new_line)
    else:
        sys.stdout.write(line)

sys.exit()
