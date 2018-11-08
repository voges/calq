# CALQ

Coverage-adaptive lossy quality value compression

[![Build Status](https://travis-ci.com/voges/calq-closed.svg?token=GduwvyBbRmwxKZ5dTnUU&branch=develop)](https://travis-ci.com/voges/calq-closed)

---

## Build instructions

Clone the CALQ repository with either

    git clone https://github.com/voges/calq.git

Build the executable from the command line with the following commands; alternatively use the CMake GUI.

    mkdir build
    cd build
    cmake ..
    make

This generates a CALQ executable named ``calq`` in the ``build`` folder.

## Usage examples

As usual, a list of the available command line options can be obtained via ``calq --help`` or ``calq -h``.

### Compression

The CALQ encoder accepts input files in the SAM format.

Basically, the following command can be used to compress the quality values from the SAM file ``file.sam``.

    calq file.sam

By default, the compressed quality values are written to the file ``file.sam.cq``. Furthermore, the CALQ encoder uses the following parameters by default:

* ``-b 10000`` | ``--blockSize 10000`` (block size in number of SAM records, i.e., alignments),
* ``-p 2`` | ``--polyploidy 2`` (sequence reads from a diploid organism are assumed),
* ``-q Illumina-1.8+`` | ``--qualityValueType Illumina-1.8+`` (quality values in the Illumina 1.8+ format (Phred+33, i.e., [0, 41] + 33) are assumed).

Thus, the above command is equivalent to the following command.

    calq -q Illumina 1.8+ -p 2 -b 10000 file.sam -o file.sam.cq

### Decompression

To perform the decompression of the file ``file.sam.cq``, the CALQ decoder requires the alignment information, namely the mapping positions (POS), the CIGAR strings, and the reference sequence name(s) (RNAME). This information can be passed to the CALQ decoder with the argument ``-s file.sam``. The switch ``-d`` invokes the decoder.

    calq -d -s file.sam file.sam.cq

By default, the reconstructed quality values are written to the file ``file.sam.cq.qual``.

Thus, the above command is equivalent to the following command.

    calq -d -s file.sam file.sam.cq -o file.sam.cq.qual

Finally, a SAM file containing the reconstructed quality values can be produced with the Python script ``replace_qual_sam.py``. This and other supplementary scripts can be found in the folder ``src/ngstools``.

    replace_qual_sam.py file.sam file.sam.cq.qual 1> file.sam.cq.sam

This produces a new SAM file ``file.sam.cq.sam`` containing the reconstructed quality values.

## Who do I talk to?

Jan Voges <[voges@tnt.uni-hannover.de](mailto:voges@tnt.uni-hannover.de)>

Mikel Hernaez <[mhernaez@illinois.edu](mailto:mhernaez@illinois.edu)>
