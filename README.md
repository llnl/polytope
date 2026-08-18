[![License: BSD](https://img.shields.io/badge/License-BSD%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Build Status](https://travis-ci.org/pbtoast/polytope.svg?branch=master)](https://travis-ci.org/pbtoast/polytope)

# Polytope

Copyright (c) 2026, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory
Written by Mike Owen, David Starinshak, and Jeffrey Johnson.
Updated by Landon Owen.
LLNL-CODE-647432
All rights reserved.

Polytope is a C++ library for generating polygonal and polyhedral meshes.
It makes use of various 2D and 3D tessellation techniques, but provides a
single representation for these tessellations, and a simple interface for
generating them.

It includes bindings for Python. These bindings allow you to easily incorporate
Polytope into your own mesh generation tools.

## Installation

Polytope works on most Linux and Mac systems.

### Software Requirements

+ A C++17 compiler
+ MPI for parallelism
+ The Bourne Again SHell (bash), for bootstrapping
+ CMake 3.21+, for configuring and generating build files
+ GNU Make or Ninja, for performing the actual build

If you want to build Python bindings, you also need the following:
+ A [Python 3](https://www.python.org/downloads) interpreter
+ [pybind11](https://github.com/pybind/pybind11), a Python/C++11
  interoperability layer
+ [PYB11Generator](https://github.com/jmikeowen/PYB11Generator), a code
  generator that processes binding definitions in Python. PYB11Generator
  produces C++ code that uses pybind11 to expose your C++ classes as Python
  classes.

### Building

To build polytope on a UNIX-like system, open `bootstrap` and modify the
variables as necessary. Then run the script using

```
./bootstrap build_dir
```

where `build_dir` is the directory in which you want to build.

Then just follow the onscreen directions: you change to that build directory,
edit `config.sh` to define your build, run it with `sh config.sh`, and start
the build using your generator's build process. For the default generator
(UNIX makefiles), this is just `make`.

### Installing

To install polytope, use the install command for the generator you've selected.
For example, if you're using a generator that writes UNIX makefiles, run

```
make install [-j #threads]
```

from your build directory.

## Python Interface

The build system generates a python virtual environment in both the build and install.
The virtual environment is in the directory PolytopePy and can be activated by

```
source PolytopePy/bin/activate
```
or
```
source PolytopePy/bin/activate.csh
```
depending on your terminal shell.

## Other Considerations

Polytope provides interfaces for a number of geometry-related tools:

+ [Triangle](http://www.cs.cmu.edu/~quake/triangle.html) by Jonathan Shewchuk
  at Berkeley
+ [Tetgen](http://www.wias-berlin.de/software/index.jsp?id=TetGen&lang=1) by
  Hang Si at Weierstrass Institute for Applied Analysis and Stochastics
+ [Boost.Polygon.Voronoi](https://www.boost.org/doc/libs/1_61_0/libs/polygon/doc/voronoi_main.htm),
  part of the Boost C++ Library

To use these tools, adjust the related CMake variables in
the config.sh that bootstrap creates. Specifically, `POLYTOPE_ENABLE_{TOOLNAME}=ON`
and `{toolname}_DIR=/path/to/tool/install`.

To use Triangle, set `POLYTOPE_ENABLE_TRIANGLE=ON` and set
`triangle_SRC_DIR` to point to the **SOURCE** of a local Triangle repo.
This is different than other TPLs where you must point to an existing TPL
install.

### Using Triangle and Tetgen

If using Triangle or Tetgen, please note: **you must comply with the licenses for these tools**.
Briefly, this means that if you want to use Triangle or Tetgen in a commercial
application, you must contact the author for permission.

To keep things simple, we don't distribute the source for either of these
tools. Only use these tools after you've made arrangements to comply with the license(s).
