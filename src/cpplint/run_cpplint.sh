#!/bin/bash
for f in $(find ../calq/ -type f -not -path "../calq/tclap/*" -not -path "../calq/Compressors/range/*"); do python2 cpplint.py --root=src --linelength=160 --extensions=c,cc,cpp,h,hh,hpp $f; done
