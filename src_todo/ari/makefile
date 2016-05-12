#
# The MIT License (MIT)
#
# Copyright (c) 2014 Mark Thomas Nelson
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# This code was written to illustrate the article:
#  Data Compression With Arithmetic Coding
#  by Mark Nelson
# published at: http://marknelson.us/2014/10/19/data-compression-with-arithmetic-coding
#
#
#CC=clang++ -std=c++1y -DLOG
CC=g++ --std=c++0x
HEADERS=compressor.h decompressor.h modelA.h model_metrics.h bitio.h byteio.h

all: tester compress decompress fp_proto

fp_proto:	fp_proto.o
	$(CC) -o fp_proto fp_proto.o

fp_proto.o:	fp_proto.cpp
	$(CC) -c fp_proto.cpp

tester: tester.o
	$(CC) -o tester tester.o

tester.o:	tester.cpp $(HEADERS)
	$(CC) -c tester.cpp

compress: compress.o
	$(CC) -o compress compress.o

compress.o:	compress.cpp $(HEADERS)
	$(CC) -DLOG -c compress.cpp

decompress: decompress.o
	$(CC) -o decompress decompress.o

decompress.o:	decompress.cpp $(HEADERS)
	$(CC) -DLOG -c decompress.cpp

clean:
	rm -f tester tester.o compress compress.o decompress.o decompress fp_proto fp_proto.o

test:
	find corpus -type f -print0 | xargs -0r -n 1 ./tester
 
