#!/usr/bin/env bash

set -eux

mkdir -p web-example
cp emscripten/web-template/* web-example
cp src/CoCoA-5/CoCoAInterpreter.wasm web-example
cp src/CoCoA-5/CoCoAInterpreter web-example/CoCoAInterpreter.js

cp src/CoCoA-5/cocoa-fs.json web-example/
cp -r src/CoCoA-5/packages web-example/
cp -r src/CoCoA-5/tests web-example/
cp -r src/CoCoA-5/CoCoAManual web-example/

cd web-example
python3 -m http.server 9999
