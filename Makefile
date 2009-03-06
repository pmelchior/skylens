INCLPATH = ./include
SRCPATH = ./src
LIBPATH = ./lib/$(SUBDIR)
LIBNAME = skylens

SHAPELENSPATH = $(HOME)/include/shapelens
NUMLAPATH = $(HOME)/include/numla
ATLASPATH = $(HOME)/include/atlas
LIBASTROPATH = $(HOME)/include/libastro

SRC = $(wildcard $(SRCPATH)/*.cc)
OBJECTS = $(SRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)
HEADERS = $(wildcard $(INCLPATH)/*.h)

CC = g++
CFLAGS = -ansi -g -Wno-deprecated -O3 -march=pentium4 -I$(INCLPATH) -I$(HOME)/include -I$(SHAPELENSPATH) -I$(NUMLAPATH) -I$(ATLASPATH) -I$(LIBASTROPATH) -DDATAPATH=$(PWD)/data

AR = ar
ARFLAGS = -sr
SHAREDFLAGS = -shared -fPIC 

.DEFAULT: all

.PHONY: clean

all: lib shared docs

clean:
	rm -f $(LIBPATH)/*

cleanlib:
	rm -f $(LIBPATH)/*

cleanshared:
	rm -f $(LIBPATH)/lib$(LIBNAME).so

lib: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).so

docs: $(HEADERS)
	doxygen Doxyfile

$(LIBPATH)/lib$(LIBNAME).a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

$(LIBPATH)/lib$(LIBNAME).so: $(OBJECTS)
	$(CC) $(SHAREDFLAGS) -o $@ $?

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@
