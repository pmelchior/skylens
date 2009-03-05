# PATHs are filled by local.py

INCLPATH = ./include
SRCPATH = ./src
LIBPATH = ./lib
PROGSRCPATH = ./progs
PROGPATH = ./bin
LIBNAME = skylens

SRC = $(wildcard $(SRCPATH)/*.cc)
OBJECTS = $(SRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)
HEADERS = $(wildcard $(INCLPATH)/*.h)
PROGS = $(wildcard $(PROGSRCPATH)/*.cc)
PROGSOBJECTS =  $(PROGS:$(PROGSRCPATH)/%.cc=$(PROGPATH)/%)

CC = g++-4.2
CFLAGS = -ansi -g -Wno-deprecated -O3 -march=pentium4 -I$(INCLPATH) -I$(HOME)/include -DDATAPATH=$(PWD)/data
CFLAG_PROGS = -L$(LIBPATH) -L$(HOME)/lib
LIBS = -l$(LIBNAME) -lastrocpp -lgsl -lgslcblas -lCCfits -lcfitsio -lfftw3

AR = ar
ARFLAGS = -sr
SHAREDFLAGS = -shared -fPIC 

.DEFAULT: all

.PHONY: clean

all: lib shared docs progs

clean:
	rm -f $(LIBPATH)/*
	rm -f $(PROGSOBJECTS)

cleanlib:
	rm -f $(LIBPATH)/*

cleanshared:
	rm -f $(LIBPATH)/lib$(LIBNAME).so

cleanprogs:
	rm -f $(PROGSOBJECTS)

lib: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).so

progs: $(PROGSOBJECTS)

docs: $(HEADERS)
	doxygen Doxyfile

$(LIBPATH)/lib$(LIBNAME).a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

$(LIBPATH)/lib$(LIBNAME).so: $(OBJECTS)
	$(CC) $(SHAREDFLAGS) -o $@ $?

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CC) $(CFLAGS) $(CFLAG_PROGS) $< -o $@ $(LIBS)

