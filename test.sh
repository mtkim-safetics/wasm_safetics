#!/bin/bash
emcc -v \
    -I src/external/Eigen \
    src/safetics/eigen.cpp \
    -lembind \
    -o wasm/eigen-sfd.js \
    --no-entry -s EXPORT_ES6=1 -s WASM=1 -s ENVIRONMENT='web' -s USE_ES6_IMPORT_META=0 -s ASSERTIONS -s MODULARIZE=1 -O1 \
    -std=c++17 \
    -sNO_DISABLE_EXCEPTION_CATCHING \
    -Wno-enum-constexpr-conversion
sed -i '' '1s;^;\/* eslint-disable *\/;' wasm/eigen-sfd.js
sed -i '' 's|eigen.wasm|/wasm/eigen.wasm|' wasm/eigen-sfd.js
perl -i -p0e "s/(if \(typeof exports === 'object' && typeof module === 'object'\))[\s\S]*/export default Module;/g" wasm/eigen-sfd.js
