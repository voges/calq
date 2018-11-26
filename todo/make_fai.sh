#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 1 ]; then printf "Usage: $0 file.[fa|fasta]\n"; exit -1; fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .fa or .fasta

printf "Input FASTA file: $1\n"

if [ ! -f $1 ]; then printf "Error: Input FASTA file $1 is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

install_path="/project/dna/install"
samtools="$install_path/samtools-1.3/bin/samtools"

if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi

###############################################################################
#                                    Make                                     #
###############################################################################

printf "Constructing FASTA index file: $1.fai\n"
if [ -f $1.fai ]; then
    printf "$1.fai already exists (not reproducing it)\n"
else
    printf "Handing over to Samtools\n"
    $samtools faidx $1
    printf "Returned from Samtools\n"
fi

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"

