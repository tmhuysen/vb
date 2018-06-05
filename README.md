# vb v0.0.1
[![Build Status](https://travis-vb.org/GQCG/vb.svg?branch=master)](https://travis-vb.org/GQCG/vb)

A C++ library for performing valence bond (vb) calculations.

## Dependencies
[![Boost Dependency](https://img.shields.io/badge/Boost-1.65.1+-000000.svg)](http://www.boost.org)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3.3.4+-000000.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-000000.svg)](https://github.com/evaleev/libint)
[![cpputil Dependency](https://img.shields.io/badge/cpputil-1.3.0+-blue.svg)](https://github.com/GQCG/cpputil)
[![bmqc Dependency](https://img.shields.io/badge/bmqc-1.2.0+-blue.svg)](https://github.com/GQCG/bmqc)
[![libwint Dependency](https://img.shields.io/badge/libwint-3.0.0+-blue.svg)](https://github.com/GQCG/libwint)
[![hf Dependency](https://img.shields.io/badge/hf-3.0.0+-blue.svg)](https://github.com/GQCG/hf)
[![numopt Dependency](https://img.shields.io/badge/numopt-1.1.0+-blue.svg)](https://github.com/GQCG/numopt)


## Installation
To install this library:
1. clone the master branch, which contains the latest release

        https://github.com/GQCG/vb.git --branch master --single-branch
        cd vb

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libvb.a` will be installed in `prefix/vb/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/vb/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(vb 0.0.1)

CMake then provides the commands `vb_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `vb` to be used in your `target_link_libraries`.
