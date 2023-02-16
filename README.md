
# C++ Quad Wavelet Tree

Wavelet trees[[1](#bib)] are a self-indexing rank and select data structure, i.e., they can answer *rank* (how often does a symbol occur in a prefix of length $i$) and *select* (where does a symbol occur for the $i$-th time) queries, while still allowing to *access* the text.

This makes them an important building block for compressed full-text indices, where they are used to answer rank queries during the pattern matching algorithm, the so called backwards search.

Using a wavelet tree requires $n\lceil\log \sigma \rceil (1+o(1))$ bits of space for a text of length $n$ over an alphabet of size $\sigma$. Access, rank, and select queries can be answered in time $O(\log\sigma)$ time.

This repository provides a fast implementation of **wavelet matrices** in C++, since they are in practice faster than wavelet trees. A companion **Rust implementation** is available [here](https://github.com/rossanoventurini/qwt/).

To improve query performance, our implementation uses a  4-ary tree instead of a binary tree as basis of the wavelet tree. The 4-ary tree layout of a wavelet tree helps to halve the number of cache misses during queries and thus reduces the query latency. Our experimental evaluation shows that our 4-ary wavelet tree can improve the latency of access, rank and select queries by factor of $\approx$ 2 compared to other implementations of wavelet trees contained in the widely used Succinct Data Structure Library ([SDSL](https://github.com/simongog/sdsl-lite)). For more details, see [Benchmarks](#bench) and the paper [[2](#bib)].

## <a name="bench">Benchmarks</a>
All the experiments are performed using a single thread on a server machine with 8 Intel~i9-9900KF cores with base frequencies of 3.60 GHz running Linux 5.19.0. Each core has a dedicated L1 cache of size 32 KiB, a dedicated L2 cache of size 256 KiB, a shared L3 cache of size 16 MiB, and 64 GiB of RAM.
The code is compiled with GCC 12.2.0 using the highest optimization setting (i.e., flags *-O3 -march=native -DNDEBUG -flto*).
A more detailed experimental evaluation can be found in [[2](#bib)]. 

### Building the code
To just clone the source code, use the following.
```bash
git clone git@github.com:MatteoCeregini/quad-wavelet-tree.git
cd quad-wavelet-tree
```

To compile 
```bash
mkdir build
cd build
cmake ..
make -j
```

If you want to run the test for quad vectors, please continue with the following commands.
```bash
./test_qvector
```

If you want to run the test for quad wavelet trees, please continue with the following commands.
```bash
./create_dataset PATH_TO_TEXT PATH_TO_DATASET
./test_qwm PATH_TO_DATASET
```

`create_dataset` reads a text file, maps it to a sequence of integers, and serializes it to disk:
* `PATH_TO_TEXT` is the path to the text you want to use in the test (something like `/texts/english.txt`).
* `PATH_TO_DATASET` is the path to the serialized sequence created by `create_dataset` (for example `/texts/english.bin`).

The resulting file is then used by `create_dataset`.

## <a name="bib">Bibliography</a>
1. Roberto Grossi, Ankur Gupta, and Jeffrey Scott Vitter. *High-order entropy-compressed text indexes.* In SODA, pages 841–850. ACM/SIAM, 2003.
2. Matteo Ceregini, Florian Kurpicz, Rossano Venturini. *Faster Wavelet Trees with Quad Vectors*. Arxiv, 2023.

----
Please cite the following paper if you use this code.
```
@misc{QWT,
  author = {Matteo Ceregini, Florian Kurpicz, Rossano Venturini},
  title = {Faster Wavelet Trees with Quad Vectors},
  publisher = {arXiv},
  year = {2023},
  doi = {...},
  url = {...}
}
```
