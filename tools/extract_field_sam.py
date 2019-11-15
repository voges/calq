#!/usr/bin/env python

import argparse
import os
import sys


def main():
    # set up the command line parser
    parser = argparse.ArgumentParser(
        description='extract field from a SAM file'
    )
    parser.add_argument(
        'file_name',
        help='file name'
    )
    parser.add_argument(
        'field_id',
        help=('field ID (in range [0,11], field_ID=11 will extract all AUX '
              'fields)')
    )
    args = parser.parse_args()

    # retrieve the command line arguments
    file_name = args.file_name
    field_id = args.field_id

    # check if the user-supplied input file exists
    if not os.path.isfile(file_name):
        sys.exit("error: this is not a file: {}".format(file_name))

    # check if the user-supplied input file can be suspected to be a SAM file
    if not file_name.endswith(".sam"):
        sys.exit("error: input file name must end with '.sam'")

    # open the input file
    file = open(file_name, 'r')

    # check if field ID is valid
    field_id = int(field_id)
    if not 0 <= field_id <= 11:
        sys.exit("error: field_id must be in range [0,11]")

    # write field to stdout
    while 1:
        output = ""
        line = file.readline().rstrip('\n')
        if not line:
            break
        if not line.startswith('@'):
            fields = line.split('\t')
            if field_id != 11:  # we will treat the AUX fields separately
                output = fields[field_id] + '\n'
            else:  # extract all AUX fields
                num_fields = len(fields)
                if num_fields > 11:  # check if we have AUX fields at all
                    for i in range(11, num_fields):
                        output += fields[i]
                        if i != (num_fields - 1):
                            output += '\t'
                    output += '\n'
                else:  # no AUX fields, so just write a line break
                    output = '\n'
        sys.stdout.write(output)

    # close the input file
    file.close()


if __name__ == "__main__":
    main()
