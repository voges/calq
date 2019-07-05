# CALQ

Coverage-adaptive lossy quality value compression

[![Build Status](https://travis-ci.org/voges/calq.svg?branch=master)](https://travis-ci.org/voges/calq)

---

## Quick start on Linux

Clone the CALQ repository:

    git clone https://github.com/voges/calq.git

Build the CALQ executable using CMake:

    mkdir build
    cd build
    cmake ..
    make

This will generate the CALQ application (named ``cip``) at ``cmake-build/bin/cip``.

## Usage examples

A list of the available command line options can be obtained via ``cip --help`` or ``cip -h``.

### Compression

The CALQ encoder accepts input files in the SAM format (https://github.com/samtools/hts-specs).

The following command can be used to compress the quality values from the SAM file ``file.sam``.

    cip file.sam

The compressed quality values are written to the file ``file.sam.cip``.

### Decompression

To perform the decompression of the file ``file.sam.cip``, the CALQ decoder requires the alignment information, namely the mapping positions, the CIGAR strings, and the reference sequence name(s). This information can be passed to the CALQ decoder with the argument ``-s file.sam``. The switch ``-d`` invokes the decoder.

    cip -d -s file.sam file.sam.cip

The reconstructed quality values are written to the file ``file.sam.cip.qual``.

Finally, a SAM file containing the reconstructed quality values can be produced with the Python script ``replace_qual_sam.py``. This and other supplementary scripts can be found at https://github.com/voges/htstools.

## Continuous integration

Commits to this repository are continuously tested on **Travis CI** (https://travis-ci.org/voges/calq). Take a look at the file ``.travis.yml`` to see what is being done on Travis' (virtual) machines.

## Build system

We use **CMake** (https://cmake.org) as build system.

## Who do I talk to?

Jan Voges <[voges@tnt.uni-hannover.de](mailto:voges@tnt.uni-hannover.de)>

Mikel Hernaez <[mhernaez@illinois.edu](mailto:mhernaez@illinois.edu)>
