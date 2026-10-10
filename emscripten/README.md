Code to allow compiling CoCoA to WASM using Emscripten.

Files:

- `build.sh`: Run as `emscripten/build.sh` from a fresh copy of CoCoA.

- `web-template`: Uses 'xterm-pty' to create a "nice" interface to the Wasm CoCoA.

See `run-web-demo.sh` as an example on how to set up a working website.

Setup Instructions:

- install emscripten:
```
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install 3.1.23
./emsdk activate 3.1.23
source ./emsdk_env.sh
```
- compile CoCoA to WASM:
From root directory of this git repository, run 
```
bash emscripten/build.sh
```
Then the script will fail due to permission issue of `check-version-defines`. Don't worry, run it again
```
bash emscripten/build.sh
```
- setup website:
From root directory of this git repository, run 
```
bash emscripten/run-web-demo.sh
```
