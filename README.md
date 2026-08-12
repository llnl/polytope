[![License: BSD](https://img.shields.io/badge/License-BSD%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Build Status](https://travis-ci.org/pbtoast/polytope.svg?branch=master)](https://travis-ci.org/pbtoast/polytope)

# Polytope

Copyright (c) 2013, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory
Written by Mike Owen, David Starinshak, and Jeffrey Johnson.
LLNL-CODE-647432
All rights reserved.

Polytope is a C++ library for generating polygonal and polyhedral meshes.
It makes use of various 2D and 3D tessellation techniques, but provides a
single representation for these tessellations, and a simple interface for
generating them.

Polytope has a simple C interface for use with other languages. It also
includes bindings for Python. These bindings allow you to easily incorporate
Polytope into your own mesh generation tools.

## Installation

Polytope works on most Linux and Mac systems.

### Software Requirements

+ A C++11 compiler (and a C compiler if you want the C interface).
+ MPI for parallelism
+ The Bourne Again SHell (bash), for bootstrapping
+ CMake 3.1+, for configuring and generating build files
+ GNU Make or Ninja, for performing the actual build

If you want to build Python bindings, you also need the following:
+ A [Python 2](https://www.python.org/downloads) interpreter
+ [pybind11](https://github.com/pybind/pybind11), a Python/C++11
  interoperability layer
+ [PYB11Generator](https://github.com/jmikeowen/PYB11Generator), a code
  generator that processes binding definitions in Python. PYB11Generator
  produces C++ code that uses pybind11 to expose your C++ classes as Python
  classes.

### Building

To build polytope on a UNIX-like system, change to your `polytope` source
directory and type the following:

```
./bootstrap build_dir
```

where `build_dir` is the directory in which you want to build.

Then just follow the onscreen directions: you change to that build directory,
edit `config.sh` to define your build, run it with `sh config.sh`, and start
the build using your generator's build process. For the default generator
(UNIX makefiles), this is just `make`. For Ninja (recommended if you have it),
it's `ninja`.

### Installing

To install polytope, use the install command for the generator you've selected.
For example, if you're using a generator that writes UNIX makefiles, run

```
make install [-j #threads]
```

from your build directory.

### Other Targets

These targets all work with Make and Ninja.

+ `test` - Runs all unit tests for the library. Use `ctest -j #threads`
   instead, though, to run the tests in parallel.
+ `clean` - Removes all build assets but retains configuration options.
+ `distclean` - Performs clean and completely removes the build directory.

## Other Considerations

Polytope provides interfaces for a number of geometry-related tools:

+ [Triangle](http://www.cs.cmu.edu/~quake/triangle.html) by Jonathan Shewchuk
  at Berkeley
+ [Tetgen](http://www.wias-berlin.de/software/index.jsp?id=TetGen&lang=1) by
  Hang Si at Weierstrass Institute for Applied Analysis and Stochastics
+ [Boost.Polygon.Voronoi](https://www.boost.org/doc/libs/1_61_0/libs/polygon/doc/voronoi_main.htm),
  part of the Boost C++ Library

### Using Triangle and Tetgen

It's easy to use Triangle and Tetgen to generate tessellations. There's only
one complication: **you must comply with the licenses for these tools**.
Briefly, this means that if you want to use Triangle or Tetgen in a commercial
application, you must contact the author for permission.

To keep things simple, we don't distribute the source for either of these
tools. Once you've made arrangements to comply with the license(s), you can
copy the source for the tool(s) into place:

+ For Triangle, copy `triangle.h` and `triangle.c` to the `src` directory.
+ For Tetgen, copy `tetgen.h` and `tetgen.cxx` to the `src` directory.

Then (re)configure and build.

