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
CFLAGS = -ansi -g $(SPECIALFLAGS) -I$(INCLPATH) -I$(NUMLAPATH) -I$(LIBASTROPATH) -I$(SHAPELENSPATH) -I$(SKYDBPATH) -DDATAPATH=$(PWD)/data
CFLAG_LIBS = -I$(HOME)/include -L$(ITALIBSLIBPATH) -L$(LIBPATH)
LIBS = -lskylens -lskydb -lshapelens -lastrocpp -lgsl -lcblas -llapack_atlas -latlas -llapack -lg2c -lCCfits -lcfitsio -lmysqlclient -lfftw3

AR = ar
ARFLAGS = -sr
SHAREDFLAGS = -shared -fPIC 

.DEFAULT: all

.PHONY: clean

all: lib shared

clean:
	rm -f $(LIBPATH)/*

cleanshared:
	rm -f $(LIBPATH)/lib$(LIBNAME).so

cleanprogs:
	rm -f $(PROGSOBJECTS)

lib: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).so

docs: $(HEADERS)
	doxygen Doxyfile

progs: $(PROGSOBJECTS)

install: lib shared
	mkdir -p $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).a $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).so $(ITALIBSLIBPATH)
	mkdir  -p $(ITALIBSPATH)/include/$(LIBNAME)
	cd $(INCLPATH) && find . -type f -name '*.h' -exec  cp --parents {} $(ITALIBSPATH)/include/$(LIBNAME)/ \; && cd ../
	mkdir -p $(PROGPATH)

$(LIBPATH)/lib$(LIBNAME).a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

$(LIBPATH)/lib$(LIBNAME).so: $(OBJECTS)
	$(CC) $(SHAREDFLAGS) -o $@ $?

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CC) $(CFLAGS) $(CFLAG_LIBS) $< -o $@ $(LIBS)

