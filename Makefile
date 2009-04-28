INCLPATH = ./include
SRCPATH = ./src
LIBPATH = ./lib/$(SUBDIR)
PROGSRCPATH = ./progs
LIBNAME = skylens

NUMLAPATH = $(ITALIBSPATH)/include/numla
LIBASTROPATH = $(ITALIBSPATH)/include/libastro
SHAPELENSPATH = $(ITALIBSPATH)/include/shapelens
SKYDBPATH = $(ITALIBSPATH)/include/skydb

SRC = $(wildcard $(SRCPATH)/*.cc)
OBJECTS = $(SRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)
HEADERS = $(wildcard $(INCLPATH)/*.h)

ITALIBSLIBPATH = $(ITALIBSPATH)/lib/$(SUBDIR)
PROGPATH = $(ITALIBSPATH)/bin/$(SUBDIR)
PROGS = $(wildcard $(PROGSRCPATH)/*.cc)
PROGSOBJECTS = $(PROGS:$(PROGSRCPATH)/%.cc=$(PROGPATH)/%)

CC = g++
CFLAGS = -ansi -bind_at_load -g $(SPECIALFLAGS) -I$(HOME)/include -I$(INCLPATH) -I$(NUMLAPATH) -I$(LIBASTROPATH) -I$(SHAPELENSPATH) -I$(SKYDBPATH) -DDATAPATH=$(PWD)/data -DHAS_FFTW3 -DSHAPELETDB=MySQL
CFLAG_LIBS = -I$(HOME)/include -L$(ITALIBSLIBPATH) -L$(LIBPATH) -L$(HOME)/lib
LIBS = -lskylens -lshapelens -lastrocpp -lgsl -lcblas -llapack_atlas -latlas -llapack -lg2c -lCCfits -lcfitsio -lmysqlclient -lfftw3 -lspatialindex

AR = ar
ARFLAGS = -sr
SHAREDFLAGS = -dynamiclib -fPIC

.DEFAULT: all

.PHONY: clean

all: lib shared

clean:
	rm -f $(LIBPATH)/*

cleanshared:
	rm -f $(LIBPATH)/lib$(LIBNAME).dylib

cleanprogs:
	rm -f $(PROGSOBJECTS)

lib: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).dylib

docs: $(HEADERS)
	doxygen Doxyfile

progs: $(PROGSOBJECTS)

install: lib shared
	mkdir -p $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).a $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).dylib $(ITALIBSLIBPATH)
	mkdir  -p $(ITALIBSPATH)/include/$(LIBNAME)
	cd $(INCLPATH) && find . -type f -name '*.h' -exec  cp --parents {} $(ITALIBSPATH)/include/$(LIBNAME)/ \; && cd ../
	mkdir -p $(PROGPATH)

$(LIBPATH)/lib$(LIBNAME).a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

$(LIBPATH)/lib$(LIBNAME).dylib: $(OBJECTS)
	$(CC) $(SHAREDFLAGS) $(CFLAG_LIBS) -o $@ $? $(LIBS)

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CC) $(CFLAGS) $(CFLAG_LIBS) $< -o $@ $(LIBS)

