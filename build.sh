#!/bin/bash

# Find all .cpp files in the directory tree
for file in $(find . -name "*.cpp"); do
  # Extract the base name of the file for the output file name
  base_name=$(basename $file .cpp)
  echo "compile $base_name.cpp to $base_name.js"
  # Create a Makefile for the file
  em++ --no-entry src/${base_name}.cpp -o src/${base_name}.js -s EXPORT_ES6=1 -s WASM=1 -s ENVIRONMENT='web' -s EXPORT_NAME='createModule' -s USE_ES6_IMPORT_META=0 -s 'EXPORTED_RUNTIME_METHODS=["ccall"]' -s ASSERTIONS -O3
done

# Find all .cpp files in the directory tree
for file in $(find . -name "*.c"); do
  # Extract the base name of the file for the output file name
  base_name=$(basename $file .c)
  echo "Building $base_name"
  # Create a Makefile for the file
  emcc --no-entry src/${base_name}.c -o src/${base_name}.js -s EXPORT_ES6=1 -s WASM=1 -s ENVIRONMENT='web' -s EXPORT_NAME='createModule' -s USE_ES6_IMPORT_META=0 -s 'EXPORTED_RUNTIME_METHODS=["ccall"]' -O3
done
