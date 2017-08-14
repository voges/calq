#!/bin/bash
for f in $(find . -type f -not -path "./calq/tclap/*" -not -path "./calq/Compressors/range/*"); do python cpplint.py --root=src --linelength=160 --extensions=c,cc,cpp,h,hh,hpp $f; done
