#!/bin/bash
git clone https://github.com/emscripten-core/emsdk.git
echo "clone emsdk done"
cd emsdk && ./emsdk install latest && ./emsdk activate latest
echo "install emsdk done"
source ./emsdk_env.sh
echo "source emsdk_env.sh done"
cd ..
