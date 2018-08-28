#!/usr/bin/env python

import sys

# Usage
if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.sam 1>clipped.sam".format(sys.argv[0]))

# Get SAM file
sam_file_name = sys.argv[1]
if not sam_file_name.endswith(".sam"):
    sys.exit("SAM file name must end with '.sam'")
sam_file = open(sam_file_name, 'r')
#print("SAM file: {}".format(sam_file_name))

sys.stderr.write("Clipping quality values into Illumina 1.8+ format ([0,41]+33)\n")

# Initialize statistics
qual_total_cnt = 0
qual_raise_cnt = 0
qual_decrease_cnt = 0

# Parse SAM file and clip the quality scores
while 1:
    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        a = line.split('\t')
        new_line = a[0] + '\t' + a[1] + '\t' + a[2] + '\t' + a[3] + '\t' + a[4] + '\t' + a[5] + '\t'
        new_line += a[6] + '\t' + a[7] + '\t' + a[8] + '\t' + a[9] + '\t'
        qual = a[10]

        # Check and clip quality scores

        # Typical ranges:
        #  Sanger         Phred+33   [0,40]
        #  Solexa         Solexa+64  [-5,40]
        #  Illumina 1.3+  Phred+64   [0,40]
        #  Illumina 1.5+  Phred+64   [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
        #  Illumina 1.8+  Phred+33   [0,41]

        # Currently using Illumina 1.8+
        new_qual = ""
        for q in qual:
            qual_total_cnt += 1
            if ord(q) < 33:
                q = chr(33)
                qual_raise_cnt += 1
            if ord(q) > 33+41:
                q = chr(33+41)
                qual_decrease_cnt += 1
            new_qual += q

        new_line += new_qual

        if len(a) > 11:
            for i in range(11, len(a)):
                new_line += '\t' + a[i]
        else:
            new_line += '\n'

        sys.stdout.write(new_line)

    else:
        sys.stdout.write(line)

sys.stderr.write("Clipped {} out of {} quality values\n".format(qual_raise_cnt+qual_decrease_cnt, qual_total_cnt))
sys.stderr.write("  Raised: {}\n".format(qual_raise_cnt))
sys.stderr.write("  Decreased: {}\n".format(qual_decrease_cnt))

sys.exit()
