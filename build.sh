#!/bin/bash

cd build && emcmake cmake .. && emmake make

# for file in *.a
# do
#   emcc -v -O3 -s WASM=1 -s SIDE_MODULE=1 -o "${file%.*}.wasm" "$file"
#   rm "$file"
# done