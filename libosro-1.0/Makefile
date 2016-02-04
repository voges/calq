CC = gcc
CFLAGS = -fPIC -O2 -Wall -Wextra -Wfloat-equal -Wshadow -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wwrite-strings -Waggregate-return -Wcast-qual -Wswitch-default -Wswitch-enum -Wunreachable-code -Wconversion -Wundef
LDFLAGS = -shared

CP = cp
MV = mv
RM = rm -f
MKDIR = mkdir -p

prefix ?= /usr/local
incdir = $(prefix)/include
libdir = $(prefix)/lib
srcdir = ./src

sources = $(shell find $(srcdir) -name '*.c')
objects = $(sources:.c=.o)

libosro = libosro.so

.PHONY: default all install clean

default: all

all: $(libosro)

install:
	$(MKDIR) $(incdir)
	$(CP) $(srcdir)/*.h $(incdir)
	$(MKDIR) $(libdir)
	$(MV) $(libosro) $(libdir)

$(libosro): $(objects)
	$(CC) ${LDFLAGS} -o $@ $^

%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	-$(RM) $(objects) $(libosro)

