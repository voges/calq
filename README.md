# README #

## To-do's ##

* Print statistics not to stderr but to a separate logfile
* Add command line argument to be able to modify the number of quantizers
* Print compression ratio for mapped as well as unmapped reads
* Add options to VariantRecalibrator: mg 4 (default is 9); -minNumBad 5000 (default is 1000)
* Option --tranche is MPEG's theta (90, 99, 99.9, 100)
* Try maximum mean discrepancy instead  of 1st- minus  2nd-largest likelihood

## What is this repository for? ##

* Adaptive lossy compression of next-generation sequencing quality values
* Current version: 1.0.0

## Build instructions

We provide a CMakeLists.txt to build the program with CMake.

### Linux + GCC ###

    mkdir build/linux_gcc
    cd build/linux_gcc
    cmake ../..
    make

### Windows + MSVC ###

    TBD

### Apple + AppleClang ###

    TBD

## Contribution guidelines ##

* Please write code conforming to the used coding style
* Please \#include with paths relative to src/. Place corresponding header
  and source files in the same directory. The source files shall
  include the corresponding header file without any path.

## Who do I talk to? ##

Jan Voges <[voges@tnt.uni-hannover.de](mailto:voges@tnt.uni-hannover.de)>
Mikel Hernaez <[mhernaez@stanford.edu](mailto:mhernaez@stanford.edu)>
