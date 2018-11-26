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
java="/usr/bin/java"
java_opts=""
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"

if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -x $picard_jar ]; then printf "Error: JAR file $picard_jar is not executable.\n"; exit -1; fi

###############################################################################
#                                    Make                                     #
###############################################################################

printf "Constructing FASTA dict file: $root.dict\n"
if [ -f $root.dict ]; then
    printf "$root.dict already exists (not reproducing it)\n"
else
    printf "Handing over to Picard\n"
    $java $java_opts -jar $picard_jar CreateSequenceDictionary R=$1 O=$root.dict
    printf "Returned from Picard\n"
fi

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"
