Clustered Elias-Fano Indexes
============================
-----------------------------
This is the code used for the experiments in the paper [*Clustered Elias-Fano Indexes*](http://pages.di.unipi.it/pibiri/papers/TOIS17.pdf), by Giulio Ermanno Pibiri and Rossano Venturini, published in ACM TOIS 2017 [1].

This guide is meant to provide a brief overview of the library and to illustrate its functionalities through some examples.
##### Table of contents
* [Building the code](#building-the-code)
* [Input data format](#input-data-format)
* [Computing the clusters](#computing-the-clusters)
* [Building the indexes](#building-the-indexes)
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


Input data format
-----------------
-----------------
The collection containing the docID and frequency lists follow the format of [`ds2i`](https://github.com/ot/ds2i), that is all integer lists are prefixed by their length written as 32-bit little-endian unsigned integers:

* `<basename>.docs` starts with a singleton binary sequence where its only
  integer is the number of documents in the collection. It is then followed by
  one binary sequence for each posting list, in order of term-ids. Each posting
  list contains the sequence of docIDs containing the term.

* `<basename>.freqs` is composed of a one binary sequence per posting list, where
  each sequence contains the occurrence counts of the postings, aligned with the
  previous file (note however that this file does not have an additional
  singleton list at its beginning).

The folder `test_data` constains an example of such collection organization. It consists in a sample of 244 postings lists drawn from Gov2 (one of the two datasets used for the experiments in the paper). For convenience all datasets have been compressed with `gzip` and must be uncompressed before running the experiments.
In particular, the `.docs` sequences have been split into two parts: these must be uncompressed and concatenated one after the other by doing

    $ cat test_collection.bin.docs.part_1 test_collection.bin.docs.part_2 \
        > test_collection.bin.docs

The folder also contains the postings lists' positions `test_collection.lists_positions.gz` and an examplar clustering `test_collection.clusters.gz` (see section [Computing the clusters](#computing-the-clusters)); a set of 500 queries named `queries`.

For the following examples, we assume to work with the sample data contained in `test_data`.

Computing the clusters
----------------------
----------------------
The executable `compute_clusters` can be used to cluster a set of postings lists, referenced from the input collection by the file listing their positions. For the other parameters of the executable, see `compute_clusters.cpp`.

As an example, the following command computes the clusters over the test collection as the ones in `test_collection.clusters.gz`:

    $ ./compute_clusters ../test_data/test_collection.bin \
                         ../test_data/test_collection.plists_positions.gz \
                         24622344 244 3 5 5 8 10 > test_collection.clusters

The computed clusters is a file listing one cluster per row. A cluster is an integer sequence: the first integer represents the number of postings lists in the cluster, the others represent the positions of the sequences belonging to the cluster. The file must be compressed with `gzip` to be used in the experiments.

Building the indexes
--------------------
--------------------
The executables `create_clustered_freq_index_fb` (frequency-based) and `create_clustered_freq_index_sb` (space-based) can be used to build clustered Elias-Fano indexes, given an input collection and a set of clusters.
For the other parameters of the executables, see the corresponding `.cpp` files. Below we show some examples.

##### Example 1.
The command

    $ ./create_clustered_freq_index_fb ../test_data/test_collection.bin \
    ../test_data/test_collection.clusters.gz 800000 clustered_opt_index.800K.bin

builds a clustered Elias-Fano index:
* using the frequency-based approach;
* whose reference list size is 800,000;
* that is serialized to the binary file `clustered_opt_index.800K.bin`.

##### Example 2.
The command

    $ ./create_freq_index opt ../test_data/test_collection.bin \
    --clusters ../test_data/test_collection.clusters.gz opt_index.bin

builds a partitioned Elias-Fano index on the same postings lists used by the corresponding clustered index (see Example 1.), as specified with the option `--clusters` and serialized to the binary file `opt_index.bin`.

##### Example 3.
The command

    $ ./create_freq_index block_interpolative ../test_data/test_collection.bin \
    --clusters ../test_data/test_collection.clusters.gz bic_index.bin

builds a Binary Interpolative index on the same postings lists used by the corresponding clustered index (see Example 1.), as specified with the option `--clusters` and serialized to the binary file `bic_index.bin`.


A comparison between the space of such indexes is summarized by the following table, where CPEF indicates the clustered Elias-Fano index, PEF the partitioned Elias-Fano index and BIC the Binary Interpolative one.

|     **Index**     |**bits x posting** |
|-------------------|-------------------|
|CPEF               |4.23               |
|PEF                |5.15 (**+17.86%**) |
|BIC                |4.60 (**+8.04%**)  |

Authors
-------
-------
* [Giulio Ermanno Pibiri](http://pages.di.unipi.it/pibiri/), <giulio.pibiri@di.unipi.it>
* [Rossano Venturini](http://pages.di.unipi.it/rossano/), <rossano.venturini@unipi.it>

Bibliography
------------
------------
* [1] Giulio Ermanno Pibiri and Rossano Venturini, *Clustered Elias-Fano Indexes*. ACM Transactions on Information Systems (TOIS 2017).
