# CALQ

Coverage-adaptive lossy quality value compression

[![build](https://travis-ci.com/voges/calq.svg?branch=master)](https://travis-ci.com/voges/calq)
[![codecov](https://codecov.io/gh/voges/calq/branch/master/graph/badge.svg)](https://codecov.io/gh/voges/calq)
[![doc](https://img.shields.io/badge/doc-online-blue)](https://voges.github.io/calq)

---

## Quick start

Clone the repository:

    git clone https://github.com/voges/calq.git

Build all libraries and executables using CMake:

    mkdir build
    cd build
    cmake ..
    make

This will generate the CALQ encoder/decoder (named ``calq-codec``) at ``build/bin/calq-codec``.

## Usage examples

A list of the available command line options can be obtained via ``calq-codec --help``.

### Compression

The CALQ encoder accepts input files in the SAM format (https://github.com/samtools/hts-specs).

The following command can be used to compress the quality values from the SAM file ``file.sam``.

    calq-codec file.sam

The compressed quality values are written to the file ``file.sam.cq``.

### Decompression

To perform the decompression of the file ``file.sam.cq`` the CALQ decoder requires the alignment information, namely the mapping positions, the CIGAR strings, and the reference sequence name(s). This information can be passed to the CALQ decoder with the argument ``-s file.sam``. The CALQ decoder uses only the alignment information from the file ``file.sam``. The switch ``-d`` invokes the decoder.

    calq-codec -d -s file.sam file.sam.cq

The reconstructed quality values are written to the file ``file.sam.cq.qual``.

Finally, a SAM file containing the reconstructed quality values can be produced with the Python script ``replace_qual_sam.py`` which is located in the ``util`` folder.

## Who do I talk to?

Jan Voges <[voges@tnt.uni-hannover.de](mailto:voges@tnt.uni-hannover.de)>

Mikel Hernaez <[mhernaez@illinois.edu](mailto:mhernaez@illinois.edu)>
