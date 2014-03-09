skylens
=========

The C++ library implements a Layer-based raytracing framework particulary well-suited for realistic simulations of weak and strong gravitational lensing. Source galaxies can be drawn from analytic models or deep spaced-based imaging. Lens planes can be populated with arbitrary deflectors, typically either from N-body simulations or analytic lens models. Both sources and lenses can be placed at freely configurable positions into the light cone, in effect allowing for multiple source and lens planes.

The library is based on a Fortran implementation by Massimo Meneghetti and has been used in several research publications, starting with [Meneghetti et al. (2008)](http://adsabs.harvard.edu/abs/2008A&A...482..403M). The principle of the C++ code is described in Section 7 of [Peter Melchior's PhD thesis](http://nbn-resolving.de/urn/resolver.pl?urn=urn:nbn:de:bsz:16-opus-109546).

If you use this library for your work, please do not forget to cite the relevant papers.

Installation
------------

The library should work under Linux and Mac OS X, with any reasonably recent C++ compiler.

### Prerequisites
* [shapelens](https://github.com/pmelchior/shapelens) as a general image handling framework. Please follow the installation instruction there. You'll need [Boost](http://www.boost.org/), [tmv](http://code.google.com/p/tmv-cpp/), and [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/).
* [GSL](http://www.gnu.org/software/gsl/)
* [SQLite3](http://www.sqlite.org)
* [FFTW](http://www.fftw.org), version >= 3.0
* [SpatialIndex](http://libspatialindex.github.io/) for fast spatial look-ups

### Compilation option

Compilations is done with a Makefile, so calling `make` in the main directory should do the trick. 

The Makefile recognizes several environment variables:
* `PREFIX`: where the library and headers are to be installed (default: `/usr/`)
* `TMV_PREFIX`: where the tmv library can be found (default: `$PREFIX`)
* `SPECIAL_FLAGS`: additional compiler flags to optimize output and provide include directories, e.g.
  `-03 -m64 -march=native -fPIC`. Make sure you add the flag `-DNDEBUG` for a production system.
* `SPECIAL_LIBS`: additional linker flags and directories, e.g.
  `-L$HOME/lib`

If the WCS library should be used, the flag `-DHAS_WCSLIB` needs to be set in `SPECIAL_FLAGS`.

`make install` copies the static and shared library to `$PREFIX/lib` and the headers to `PREFIX/include/skylens`. `make progs` compiles any `.cc` file in the `progs` directory, `make installprogs` copies the executables to `$PREFIX/bin`.

`make clean` and `make cleanprogs` delete all compiled code. 

Documentation
-------------

`make docs` runs [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) on the headers, the files go to `doc/html`. An online version of the latest development branch is available [here](http://www.physics.ohio-state.edu/~melchior.12/docs/skylens/classes.html).

