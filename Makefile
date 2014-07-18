# Makefile assumes PREFIX (root path for installing files) to be set.
# If it is not set, it will point to /usr/.
# The path to TMV is taken from TMV_PREFIX, which defaults to PREFIX.
# Additional compilation flags can be specified with SPECIAL_FLAGS, 
# linker flags with SPECIAL_LIBS

LIBNAME = skylens
INCLPATH = include
SRCPATH = src
LIBPATH = lib
DOCPATH = doc
PROGSRCPATH = progs
PROGPATH = bin

########
# setting defaults and other environment variables
########
ifndef PREFIX
	PREFIX = /usr
endif

ifndef TMV_PREFIX
	TMV_PREFIX = $(PREFIX)
endif

# which OS
UNAME := $(shell uname)

# git revision macro
GITREV = $(shell git describe --always)

# compilation flags
CC = g++
CFLAGS = -ansi $(SPECIAL_FLAGS) -DGITREV=$(GITREV)
ifneq ($(UNAME),Linux)
	CFLAGS += -bind_at_load
endif

# linker flags
TMVLINK := $(shell cat ${TMV_PREFIX}/share/tmv/tmv-link)
LIBS = -lshapelens $(TMVLINK) -lgsl -lcfitsio -lsqlite3 -lfftw3 -lspatialindex

# needed because libshapelens may have been compiled with it
ifneq (,$(findstring HAS_WCSLIB,$(SPECIAL_FLAGS)))
	LIBS += -lwcs
endif

# archiver flags
AR = ar
ARFLAGS = -sr
ifeq ($(UNAME),Linux)
	SHAREDFLAGS = -shared -fPIC 
	LIBEXT = so
else
	SHAREDFLAGS = -dynamiclib -fPIC
	LIBEXT = dylib
endif

########
# target section
########
SRC = $(wildcard $(SRCPATH)/*.cc)
OBJECTS = $(SRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)
PROGS = $(wildcard $(PROGSRCPATH)/*.cc)
PROGSOBJECTS = $(PROGS:$(PROGSRCPATH)/%.cc=$(PROGPATH)/%)
HEADERS = $(wildcard $(INCLPATH)/*.h)

all: $(LIBPATH) $(DOCPATH) library shared

$(LIBPATH):
	mkdir -p $(LIBPATH)

$(PROGPATH):
	mkdir -p $(PROGPATH)

$(DOCPATH):
	mkdir -p $(DOCPATH)

.PHONY: clean

clean: 
	rm -f $(OBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)	

cleandocs:
	rm -rf $(DOCPATH)/*

cleanprogs:
	rm -f $(PROGSOBJECTS)

library: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

install: library shared
	mkdir -p $(PREFIX)/lib
	cp $(LIBPATH)/lib$(LIBNAME).a $(PREFIX)/lib/
	cp $(LIBPATH)/lib$(LIBNAME).$(LIBEXT) $(PREFIX)/lib/
	mkdir  -p $(PREFIX)/include/$(LIBNAME)
	cp $(INCLPATH)/*.h $(PREFIX)/include/$(LIBNAME)/

progs: $(PROGPATH) $(PROGSOBJECTS)

installprogs: progs
	mkdir -p $(PREFIX)/bin
	cp $(PROGSOBJECTS) $(PREFIX)/bin/

docs: $(HEADERS)
	doxygen doc/Doxyfile


$(LIBPATH)/lib$(LIBNAME).a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

$(LIBPATH)/lib$(LIBNAME).$(LIBEXT): $(OBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)	
ifeq ($(UNAME),Linux)
	$(CC) $(SHAREDFLAGS) -o $@ $^
else
	$(CC) $(SHAREDFLAGS) $(SPECIAL_LIBS) -o $@ $^ $(LIBS)
endif

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CC) $(CFLAGS) $(SPECIAL_LIBS) -L$(LIBPATH) $< -o $@ -l$(LIBNAME) $(LIBS)
