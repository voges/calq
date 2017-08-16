#!/usr/bin/env python

###############################################################################
#             Set quality values in a SAM file to a constant                  #
###############################################################################

import sys


if len(sys.argv) != 4:
    sys.exit("Usage: python {} file.sam dc_quality_value quality_value_offset".format(sys.argv[0]))

sam_file_name = sys.argv[1]
samFile = open(sam_file_name, 'r')
dc_quality_value = sys.argv[2]
quality_value_offset = sys.argv[3]
new_sam_file_name = sam_file_name + '.dc_quality_value_{}+{}.sam'.format(quality_value_offset, dc_quality_value)
new_sam_file = open(new_sam_file_name, 'w')

print("SAM file: {}".format(sam_file_name))
print("DC quality value: {}".format(dc_quality_value))
print("Quality value offset: {}".format(quality_value_offset))
print("New SAM file: {}".format(new_sam_file_name))

nr_lines = 0
nr_header_lines = 0
nr_alignment_lines = 0

print("Setting quality values")
while 1:
    print("\r{}".format(nr_lines)),
    sys.stdout.flush()

    line = samFile.readline()
    if not line:
        break
    if line[0] != '@':
        a = line.split('\t')
        new_line = a[0] + '\t' + a[1] + '\t' + a[2] + '\t' + a[3] + '\t' + a[4] + '\t' + a[5] + '\t'
        new_line += a[6] + '\t' + a[7] + '\t' + a[8] + '\t' + a[9] + '\t'

        new_quality_values = [int(dc_quality_value+quality_value_offset) for x in range(len(a[9]))]
        new_line += ''.join(chr(q) for q in new_quality_values)

        if len(a) > 11:
            for i in range(11, len(a)):
                new_line += '\t' + a[i]
        else:
            new_line += '\n'
        new_sam_file.write(new_line)
        nr_alignment_lines += 1
    else:
        new_sam_file.write(line)
        nr_header_lines += 1
    nr_lines += 1
print("")
print("Processed {} lines ({} header lines, {} alignment lines)".format(nr_lines, nr_header_lines, nr_alignment_lines))

sys.exit()

