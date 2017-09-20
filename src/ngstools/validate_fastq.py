#!/usr/bin/env python

import atexit
import sys


error_cnt = 0
error_cnt_max = 20
error_string = ""
warning_cnt = 0
warning_cnt_max = 100
warning_string = ""

record_cnt = 0
line_cnt = 0


def error(line, record, message):
    global error_string
    global error_cnt
    error_string += "Error in line {} (record {}): {}\n".format(line, record, message)
    error_cnt += 1
    if error_cnt > error_cnt_max:
        sys.exit("Too many errors")


def warning(line, record, message):
    global warning_string
    global warning_cnt
    warning_string += "Warning in line {} (record {}): {}\n".format(line, record, message)
    warning_cnt += 1
    if warning_cnt > warning_cnt_max:
        sys.exit("Too many warnings")


def print_summary():
    global error_cnt
    global error_string
    global warning_cnt
    global warning_string
    global record_cnt
    global line_cnt
    print "Found {} error(s) and {} warning(s)".format(error_cnt, warning_cnt)
    if len(error_string) != 0:
        print error_string,
    if len(warning_string) != 0:
        print warning_string,
    print "Processed {} FASTQ record(s) (i.e. {} lines)".format(record_cnt, line_cnt)


# Usage
if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.[fastq|fq]".format(sys.argv[0]))

# Get FASTQ file
fastq_file_name = sys.argv[1]
if not fastq_file_name.endswith(".fastq") and not fastq_file_name.endswith(".fq"):
    sys.exit("FASTQ file must end with either '.fastq' or '.fq'")
fastqFile = open(fastq_file_name, 'r')
print "FASTQ file: {}".format(fastq_file_name)

atexit.register(print_summary)

allowed_raw_sequence_symbols = "ACGTNacgtn"
print "Allowed raw sequence symbols: {}".format(allowed_raw_sequence_symbols)
print "Will exit after {} errors or {} warnings".format(error_cnt_max, warning_cnt_max)

while 1:
    # Try to read a whole record (i.e. 4 lines)
    sequence_identifier = fastqFile.readline().rstrip('\n')
    if not sequence_identifier:
        break  # reached end of file, everything ok
    line_cnt += 1
    sequence_identifier_cnt = line_cnt

    raw_sequence = fastqFile.readline().rstrip('\n')
    if not raw_sequence:
        error(line_cnt+1, record_cnt+1, "Record is not complete")
    line_cnt += 1
    raw_sequence_cnt = line_cnt

    description = fastqFile.readline().rstrip('\n')
    if not description:
        error(line_cnt+1, record_cnt+1, "Record is not complete")
    line_cnt += 1
    description_cnt = line_cnt

    quality_scores = fastqFile.readline().rstrip('\n')
    if not quality_scores:
        error(line_cnt+1, record_cnt+1, "Record is not complete")
    line_cnt += 1
    quality_scores_cnt = line_cnt

    record_cnt += 1

    # Check sequence identifier
    if len(sequence_identifier) > 0:
        if len(sequence_identifier) < 2:
            error(sequence_identifier_cnt, record_cnt, "Sequence identifier is too short")
        if sequence_identifier[0] != '@':
            error(sequence_identifier_cnt, record_cnt, "Sequence identifier does not begin with @")
    else:
        error(sequence_identifier_cnt, record_cnt, "Sequence identifier is empty")

    # Check raw sequence
    if len(raw_sequence) == 0:
        error(raw_sequence_cnt, record_cnt, "Raw Sequence has length 0")
    else:
        for s in raw_sequence:
            if s not in allowed_raw_sequence_symbols:
                error(raw_sequence_cnt, record_cnt, "Invalid character ('{}') in raw sequence".format(s))

    # Check description
    if len(description) > 0:
        if description[0] != '+':
            error(description_cnt, record_cnt, "Description does not begin with '+'")
        if len(description) > 1:
            if description[1:] != sequence_identifier[1:]:
                error(description_cnt, record_cnt, "Description does not match sequence identifier")
    else:
        error(description_cnt, record_cnt, "Description is empty")

    # Check quality scores

    # Typical ranges:
    #  Sanger         Phred+33   [0,40]
    #  Solexa         Solexa+64  [-5,40]
    #  Illumina 1.3+  Phred+64   [0,40]
    #  Illumina 1.5+  Phred+64   [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
    #  Illumina 1.8+  Phred+33   [0,41]

    if len(quality_scores) > 0:
        if len(quality_scores) != len(raw_sequence):
            error(quality_scores_cnt, record_cnt, "Quality score length does not match raw sequence length")
        for q in quality_scores:
            if ord(q) < 33 or ord(q) > 126:
                error(quality_scores_cnt, record_cnt, "Invalid character ('{}'={}) in quality string".format(q, ord(q)))
            if ord(q) > 104:
                warning(quality_scores_cnt, record_cnt, "Unusually high quality score ('{}'={})".format(q, ord(q)))
    else:
        error(quality_scores_cnt, record_cnt, "Quality score line is empty")

sys.exit()
