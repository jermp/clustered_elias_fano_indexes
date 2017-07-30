Clustered Elias-Fano Indexes
============================
-----------------------------
This is the code used for the experiments in the paper [*Clustered Elias-Fano Indexes*](http://pages.di.unipi.it/pibiri/papers/TOIS17.pdf), by Giulio Ermanno Pibiri and Rossano Venturini, published in ACM TOIS 2017 [1].

This guide is meant to provide a brief overview of the library and to illustrate its funtionalities through some examples.
##### Table of contents
* [Building the code](#building-the-code)
* [Authors](#authors)
* [Bibliography](#bibliography)

Building the code
-----------------
-----------------
The code is tested on Linux Ubuntu with `gcc` 5.4.1. The following dependencies are needed for the build: `CMake` >= 2.8 and `Boost` >= 1.58.

The code is largely based on the [`ds2i`](https://github.com/ot/ds2i) project, so it depends on several submodules. If you have cloned the repository without `--recursive`, you will need to perform the following commands before
building:

    $ git submodule init
    $ git submodule update

To build the code on Unix systems (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following:

    $ mkdir build
    $ cd build
    $ cmake .. -DCMAKE_BUILD_TYPE=Release
    $ make -j[number of jobs]

Setting `[number of jobs]` is recommended, e.g., `make -j4`.

Unless otherwise specified, for the rest of this guide we assume that we type the terminal commands of the following examples from the created directory `build`.

Authors
-------
-------
* [Giulio Ermanno Pibiri](http://pages.di.unipi.it/pibiri/), <giulio.pibiri@di.unipi.it>
* [Rossano Venturini](http://pages.di.unipi.it/rossano/), <rossano.venturini@unipi.it>

Bibliography
------------
------------
* [1] Giulio Ermanno Pibiri and Rossano Venturini, *Clustered Elias-Fano Indexes*. ACM Transactions on Information Systems (TOIS 2017).
