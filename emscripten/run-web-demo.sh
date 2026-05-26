#!/usr/bin/env bash

set -eux

mkdir -p web-example
cp emscripten/web-template/* web-example
cp src/CoCoA-5/CoCoAInterpreter.wasm web-example
cp src/CoCoA-5/CoCoAInterpreter.js web-example
cp src/CoCoA-5/CoCoAInterpreter.data web-example

cd web-example
python3 -m http.server 9999
