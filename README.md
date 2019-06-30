# CALQ

Coverage-adaptive lossy quality value compression

[![Build Status](https://travis-ci.org/voges/calq.svg?branch=master)](https://travis-ci.org/voges/calq)

---

## Quick start on Linux

Clone the CALQ repository with

    git clone https://github.com/voges/calq.git

Build the CALQ application:

    mkdir build
    cd build
    cmake ..
    make

This generates an executable named ``calq`` in the ``build`` folder.

## Usage examples

A list of the available command line options can be obtained via ``calq --help`` or ``calq -h``.

### Compression

The CALQ encoder accepts input files in the SAM format (https://github.com/samtools/hts-specs).

The following command can be used to compress the quality values from the SAM file ``file.sam``.

    calq file.sam

By default, the compressed quality values are written to the file ``file.sam.cq``. Furthermore, the CALQ encoder uses the following parameters by default:

* ``-b 10000`` | ``--blockSize 10000`` (block size in number of SAM records, i.e., alignments),
* ``-p 2`` | ``--polyploidy 2`` (sequence reads from a diploid organism are assumed),
* ``-q Illumina-1.8+`` | ``--qualityValueType Illumina-1.8+`` (quality values in the Illumina 1.8+ format (Phred+33, i.e., [0, 41] + 33) are assumed).

Thus, the above command is equivalent to the following command.

    calq -q Illumina 1.8+ -p 2 -b 10000 file.sam -o file.sam.cq

### Decompression

To perform the decompression of the file ``file.sam.cq``, the CALQ decoder requires the alignment information, namely the mapping positions, the CIGAR strings, and the reference sequence name(s). This information can be passed to the CALQ decoder with the argument ``-s file.sam``. The switch ``-d`` invokes the decoder.

    calq -d -s file.sam file.sam.cq

By default, the reconstructed quality values are written to the file ``file.sam.cq.qual``.

Thus, the above command is equivalent to the following command.

    calq -d -s file.sam file.sam.cq -o file.sam.cq.qual

Finally, a SAM file containing the reconstructed quality values can be produced with the Python script ``replace_qual_sam.py``. This and other supplementary scripts can be found at https://github.com/voges/htstools.

## Continuous integration

Commits to this repository are continuously tested on **Travis CI** (https://travis-ci.org/voges/calq). Take a look at ``.travis.yml`` to see what is being done on Travis' (virtual) machines.

## Build system

We use **CMake** (https://cmake.org/) as build system and we provide a ``CMakeLists.txt`` to build CALQ.

## Version control system

### Branching

We use **Git** and we use the **Gitflow** workflow (https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).

That means:

* The ``master`` branch contains only *release* commits.
* Every commit on the master branch is *tagged* according to **Semantic Versioning 2.0.0** (see below).
* Development generally takes place on the ``develop`` branch.
* Actual development takes place in *feature* branches, e.g., ``feature/my_fancy_feature``.
* Once a *feature* is completed, its branch can be merged back into the ``develop`` branch.

### Version numbers

We use the Semantic Versioning 2.0.0 (https://semver.org).

That means:

* The **release** version number format is: MAJOR.MINOR.PATCH
* Increment the
  * MAJOR version when you make incompatible API changes,
  * MINOR version when you add functionality in a backwards-compatible manner, and
  * PATCH version when you make backwards-compatible bug fixes.
* **Pre-release** versions are denoted by appending a hyphen and a series of dot separated identifiers immediately following the patch version.
  * Example 1: 1.0.0-alpha.1 ("alpha version 1 of the planned major release 1")
  * Example 2: 1.0.0-beta.1 ("beta version 1 of the planned major release 1")
  * Example 3: 1.0.0-rc.1 ("release candidate (rc) 1 of the planned major release 1")
  * Note: 1.0.0-alpha.1 < 1.0.0-beta.1 < 1.0.0-rc.1 by definition (see https://semver.org)
 
## Who do I talk to?

Jan Voges <[voges@tnt.uni-hannover.de](mailto:voges@tnt.uni-hannover.de)>

Mikel Hernaez <[mhernaez@illinois.edu](mailto:mhernaez@illinois.edu)>
