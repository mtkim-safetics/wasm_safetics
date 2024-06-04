#!/bin/bash
emcc -v -I src/external/Eigen \
src/safetics/safetics.cpp \
--no-entry -s EXPORT_ES6=1 -s WASM=1 -s ENVIRONMENT='web' -s EXPORT_NAME='safatics' -s USE_ES6_IMPORT_META=0 -s ASSERTIONS -O2 \


