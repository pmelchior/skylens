INCLPATH = include
SRCPATH = src
LIBPATH = lib
DOCPATH = doc
PROGSRCPATH = progs
LIBNAME = skylens

SRC = $(wildcard $(SRCPATH)/*.cc)
OBJECTS = $(SRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)
HEADERS = $(wildcard $(INCLPATH)/*.h)

ITALIBSLIBPATH = $(ITALIBSPATH)/lib
PROGPATH = $(ITALIBSPATH)/bin
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

# SVN revision macro
SVNREV = $(shell svnversion -n)

# compilation flags
CFLAGS = -ansi -g $(SPECIALFLAGS) -DSVNREV=$(SVNREV)

ifneq ($(UNAME),Linux)
	CFLAGS = $(CFLAGS) -bind_at_load
endif

# libraries
ifneq (,$(findstring HAS_ATLAS,$(SPECIALFLAGS)))
	LIBS = -lskylens -lshapelens -lastrocpp -lgsl -llapack_atlas -latlas -llapack -lCCfits -lcfitsio -lsqlite3 -lfftw3 -lspatialindex
else
	LIBS = -lskylens -lshapelens -lastrocpp -lgsl -lgslcblas -lCCfits -lcfitsio -lsqlite3 -lfftw3 -lspatialindex
endif

AR = ar
ARFLAGS = -sr
ifeq ($(UNAME),Linux)
SHAREDFLAGS = -shared -fPIC 
LIBEXT = so
else
SHAREDFLAGS = -dynamiclib -fPIC
LIBEXT = dylib
endif

all: $(LIBPATH) $(DOCPATH) library shared

$(LIBPATH):
	mkdir -p $(LIBPATH)

$(DOCPATH):
	mkdir -p $(DOCPATH)

.PHONY: clean

clean:
	rm -f $(LIBPATH)/*

cleanshared:
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

cleanprogs:
	rm -f $(PROGSOBJECTS)

library: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

docs: $(HEADERS)
	doxygen Doxyfile

progs: $(PROGSOBJECTS)

install: library shared
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
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)
	$(CCOMPILER) $(SHAREDFLAGS) -o $@ $^
else
$(LIBPATH)/lib$(LIBNAME).$(LIBEXT): $(OBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)
	$(CCOMPILER) $(SHAREDFLAGS) $(CFLAG_LIBS) -o $@ $^ $(LIBS)
endif

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CCOMPILER) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CCOMPILER) $(CFLAGS) $(CFLAG_LIBS) $< -o $@ $(LIBS)

