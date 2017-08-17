#!/usr/bin/env python

import sys

# Usage
if len(sys.argv) != 3:
    sys.exit("Usage: python {} file.sam field_id".format(sys.argv[0]))

# Get SAM file
sam_file_name = sys.argv[1]
if not sam_file_name.endswith(".sam"):
    sys.exit("SAM file name must end with '.sam'")
sam_file = open(sam_file_name, 'r')

# Get field ID
field_id = int(sys.argv[2])
if not 0 <= field_id <= 11:
    sys.exit("field_id must be in range [0,11]")

# Write part to stdout
while 1:
    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        field = fields[field_id]
        sys.stdout.write(field + "\n")

sys.exit()
