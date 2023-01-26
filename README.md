### tl;dr
To just clone the source code, use the following.
```bash
git clone git@github.com:MatteoCeregini/quad-wavelet-tree.git
```
If you want to run the test for quad vectors, please continue with the following commands.
```bash
cd quad-wavelet-tree
mkdir build
cd build
cmake ..
make -j test_qvector
./test_qvector
```
If you want to run the test for quad wavelet trees, please continue with the following commands.
```bash
cd quad-wavelet-tree
mkdir build
cd build
cmake ..
make -j create_dataset
make -j test_qwm
./create_dataset PATH_TO_TEXT PATH_TO_DATASET
./test_qwm PATH_TO_DATASET
```
`create_dataset` reads a text file, maps it to a sequence of integers, and serializes it to disk:
* `PATH_TO_TEXT` is the path to the text you want to use in the test (something like `/texts/english.txt`).
* `PATH_TO_DATASET` is the path to the serialized sequence created by `create_dataset`.

The resulting file is then used by `create_dataset`.
