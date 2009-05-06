INCLPATH = ./include
SRCPATH = ./src
LIBPATH = ./lib/$(SUBDIR)
PROGSRCPATH = ./progs
LIBNAME = skylens

NUMLAPATH = $(ITALIBSPATH)/include/numla
LIBASTROPATH = $(ITALIBSPATH)/include/libastro

SRC = $(wildcard $(SRCPATH)/*.cc)
OBJECTS = $(SRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)
HEADERS = $(wildcard $(INCLPATH)/*.h)

ITALIBSLIBPATH = $(ITALIBSPATH)/lib/$(SUBDIR)
PROGPATH = $(ITALIBSPATH)/bin/$(SUBDIR)
PROGS = $(wildcard $(PROGSRCPATH)/*.cc)
PROGSOBJECTS = $(PROGS:$(PROGSRCPATH)/%.cc=$(PROGPATH)/%)

# define compiler
ifdef CCOMPILER
$(warning Using compiler $(CCOMPILER))
else
CCOMPILER = g++
endif

# which OS
UNAME := $(shell uname)

# compilation flags
CFLAGS = -ansi -g $(SPECIALFLAGS) -I$(LIBASTROPATH) -DDATAPATH=$(PWD)/data -DHAS_FFTW3 -DSHAPELETDB=MySQL

ifneq ($(UNAME),Linux)
	CFLAGS = $(CFLAGS) -bind_at_load
endif

# flags for linking
#CFLAG_LIBS = -L$(ITALIBSLIBPATH) -L$(LIBPATH)
# libraries
LIBS = -lskylens -lshapelens -lastrocpp -lskydb -lgsl -lcblas -llapack_atlas -latlas -llapack -lg2c -lCCfits -lcfitsio -lmysqlclient -lfftw3 -lspatialindex

AR = ar
ARFLAGS = -sr
ifeq ($(UNAME),Linux)
SHAREDFLAGS = -shared -fPIC 
LIBEXT = so
else
SHAREDFLAGS = -dynamiclib -fPIC
LIBEXT = dylib
endif

.DEFAULT: all

.PHONY: clean

all: lib shared

clean:
	rm -f $(LIBPATH)/*

cleanshared:
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

cleanprogs:
	rm -f $(PROGSOBJECTS)

lib: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

docs: $(HEADERS)
	doxygen Doxyfile

progs: $(PROGSOBJECTS)

install: lib shared
	mkdir -p $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).a $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).$(LIBEXT) $(ITALIBSLIBPATH)
	mkdir  -p $(ITALIBSPATH)/include/$(LIBNAME)
	cd $(INCLPATH) && find . -type f -name '*.h' -exec  cp --parents {} $(ITALIBSPATH)/include/$(LIBNAME)/ \; && cd ../
	mkdir -p $(PROGPATH)

$(LIBPATH)/lib$(LIBNAME).a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

ifeq ($(UNAME),Linux)
$(LIBPATH)/lib$(LIBNAME).$(LIBEXT): $(OBJECTS)
	$(CCOMPILER) $(SHAREDFLAGS) -o $@ $?
else
$(LIBPATH)/lib$(LIBNAME).$(LIBEXT): $(OBJECTS)
	$(CCOMPILER) $(SHAREDFLAGS) $(CFLAG_LIBS) -o $@ $? $(LIBS)
endif

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CCOMPILER) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CCOMPILER) $(CFLAGS) $(CFLAG_LIBS) $< -o $@ $(LIBS)

