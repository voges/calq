#!/usr/bin/env python

###############################################################################
#          Uniform quantization of quality values in a FASTQ file             #
###############################################################################

import sys


if len(sys.argv) != 4:
    sys.exit("Usage: python {} file.fastq qmin:qmax nr_bins".format(sys.argv[0]))

# Get FASTQ file
fastq_file_name = sys.argv[1]
if not fastq_file_name.endswith(".fastq") and not fastq_file_name.endswith(".fq"):
    sys.exit("FASTQ file name must end with either '.fastq' or '.fq'")
fastq_file = open(fastq_file_name, 'r')
print("FASTQ file: {}".format(fastq_file_name))

# Get min/max quality values
minMax = sys.argv[2].split(":")
qmin = int(minMax[0])
qmax = int(minMax[1])
qrange = qmax - qmin
print("Quality value range: [{},{}]".format(qmin, qmax))

# Get number of quantization bins
nr_bins = int(sys.argv[3])
if nr_bins < 2 or nr_bins > 8:
    sys.exit("nr_bins must be in [2,8]")
print("Number of bins: {}".format(nr_bins))

# Open output file
out_file_name = fastq_file_name + ".{}-qmin{}-qmax{}-num_bins{}.fq".format(sys.argv[0], qmin, qmax, nr_bins)
out_file = open(out_file_name, 'w')
print("Out file: {}".format(out_file_name))

# Compute quantizer step size
step_size = int(qrange / nr_bins)
print("Step size: {}".format(step_size))

# Compute bounds
print("Bounds:")
bounds = {}
for b in range(0, nr_bins + 1):
    bounds[b] = qmin + b * step_size
if bounds[nr_bins] != qmax:
    #   bounds[nr_bins] += 1
    bounds[nr_bins] = qmax
for b in range(0, nr_bins + 1):
    print("  {} -> {} ".format(b, bounds[b]))

# Compute representative values
print("Representative values:")
representative_values = {}
for r in range(0, nr_bins):
    representative_values[r] = bounds[r] + int(step_size / 2)
    print("  {} -> {}".format(r, representative_values[r]))

# Compute quantization table
print("Quantization table:")
quantization_table = {}
idx = 0
for q in range(qmin, qmax + 1):
    if q > bounds[idx + 1]:
        idx += 1
    quantization_table[q] = representative_values[idx]
    print("  {} -> {}".format(q, quantization_table[q]))

record_cnt = 0
line_cnt = 0

while 1:
    sequence_identifier = fastq_file.readline().rstrip('\n')
    if not sequence_identifier:
        break  # reached end of file, everything ok
    out_file.write(sequence_identifier + '\n')
    line_cnt += 1

    raw_sequence = fastq_file.readline().rstrip('\n')
    out_file.write(raw_sequence + '\n')
    line_cnt += 1

    description = fastq_file.readline().rstrip('\n')
    out_file.write(description + '\n')
    line_cnt += 1

    quality_scores = fastq_file.readline().rstrip('\n')
    line_cnt += 1

    record_cnt += 1

    # Check quality score line
    if len(quality_scores) > 0:
        for q in quality_scores:
            qq = quantization_table[ord(q)]
            out_file.write(chr(qq))
            # print("'{}' = {} -> {}".format(q, ord(q), qq))
        out_file.write('\n')
    else:
        sys.exit("Quality score line (record {}, line {}) is empty".format(record_cnt, line_cnt))

print("Processed {} record(s) (i.e. {} lines)".format(record_cnt, line_cnt))

out_file.close()

sys.exit()

