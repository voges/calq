# CALQ

Coverage-adaptive lossy quality value compression

[![travis-ci](https://travis-ci.org/voges/calq.svg?branch=master)](https://travis-ci.org/voges/calq)
[![codecov](https://codecov.io/gh/voges/calq/branch/master/graph/badge.svg)](https://codecov.io/gh/voges/calq)
[![doc](https://img.shields.io/badge/doc-online-blue)](https://voges.github.io.calq)

---

## Quick start

Clone the repository:

    git clone https://github.com/voges/calq.git

Build all libraries and executables using CMake:

    mkdir build
    cd build
    cmake ..
    make

This will generate the CALQ encoder/decoder (named ``cip``) at ``build/bin/cip``.

## Usage examples

A list of the available command line options can be obtained via ``cip --help``.

### Compression

The CALQ encoder accepts input files in the SAM format (https://github.com/samtools/hts-specs).

The following command can be used to compress the quality values from the SAM file ``file.sam``.

    cip file.sam

The compressed quality values are written to the file ``file.sam.cip``.

### Decompression

To perform the decompression of the file ``file.sam.cip`` the CALQ decoder requires the alignment information, namely the mapping positions, the CIGAR strings, and the reference sequence name(s). This information can be passed to the CALQ decoder with the argument ``-s file.sam``. The switch ``-d`` invokes the decoder.

    cip -d -s file.sam file.sam.cip

The reconstructed quality values are written to the file ``file.sam.cip.qual``.

Finally, a SAM file containing the reconstructed quality values can be produced with the Python script ``replace_qual_sam.py`` which is located in the ``scripts`` folder.

## Development

### Continuous integration

Commits to this repository are continuously tested on **Travis CI** (https://travis-ci.org/voges/calq). Take a look at the file ``.travis.yml`` to see what is being done on Travis' (virtual) machines.

### Documentation

The Doxygen documentation of the CALQ library (i.e., the documentation of all source code in ``source/libs/calq/``) is available at https://voges.github.io/calq/.

### Code coverage

Code coverage reports are available at https://codecov.io/gh/voges/calq.

## Who do I talk to?

Jan Voges <[voges@tnt.uni-hannover.de](mailto:voges@tnt.uni-hannover.de)>

Mikel Hernaez <[mhernaez@illinois.edu](mailto:mhernaez@illinois.edu)>
