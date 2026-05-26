#!/usr/bin/env bash

set -eux

# Executed from the root of the git repo
BASEDIR="$(pwd)"
EXTERN_DIR="$BASEDIR/extern"

if ! command -v emmake &> /dev/null; then
    echo "Please install, and source, emscripten"
    echo "See https://emscripten.org/docs/getting_started/downloads.html for install instructions"
    exit 1
fi

AUX_BUILD="$EXTERN_DIR/emscripten/build"
AUX_PREFIX="$EXTERN_DIR/emscripten/install"

mkdir -p "$AUX_BUILD"
mkdir -p "$AUX_PREFIX"
mkdir -p "$AUX_PREFIX/lib"
mkdir -p "$AUX_PREFIX/include"
mkdir -p "$EXTERN_DIR"

# ------------------------------------------------------------------------------
# GMP 
# ------------------------------------------------------------------------------
(
    mkdir -p "$AUX_BUILD/gmp"
    cd "$AUX_BUILD/gmp"
    if [[ ! -d "$EXTERN_DIR/gmp" ]]; then
        echo "Downloading GMP source..."
        hg clone https://gmplib.org/repo/gmp/ "$EXTERN_DIR/gmp"
        cd "$EXTERN_DIR/gmp" && ./.bootstrap && cd -
    fi
    if [[ ! -f config.status ]]; then
        CC_FOR_BUILD=/usr/bin/gcc ABI=standard \
        emconfigure "$EXTERN_DIR/gmp/configure" \
            --build i686-pc-linux-gnu --host none \
            --disable-assembly --enable-cxx \
            --prefix="$AUX_PREFIX"
    fi
    emmake make -j8
    emmake make install
)

# ------------------------------------------------------------------------------
# BOOST
# ------------------------------------------------------------------------------
(
    BOOST_VERSION_DIR="boost_1_91_0"
    BOOST_VERSION_URL="1.91.0"
    mkdir -p "$AUX_BUILD/boost"
    cd "$AUX_BUILD/boost"
    
    if [[ ! -d "$EXTERN_DIR/$BOOST_VERSION_DIR" ]]; then
        echo "Downloading Boost $BOOST_VERSION_URL..."
        curl -sSL "https://archives.boost.io/release/$BOOST_VERSION_URL/source/${BOOST_VERSION_DIR}.tar.gz" | tar xz -C "$EXTERN_DIR"
    fi
    
    cd "$EXTERN_DIR/$BOOST_VERSION_DIR"
    
    if [[ ! -f b2 ]]; then
        env -u CC -u CXX -u CFLAGS -u CXXFLAGS -u LDFLAGS \
        ./bootstrap.sh --with-toolset=gcc
    fi
    
    cat > user-config-wasm.jam <<EOF
using emscripten : : em++ : 
    <cxxflags>"${CXXFLAGS:-}" 
    <cflags>"${CFLAGS:-}" 
    <linkflags>"${LDFLAGS:-}" 
    <archiver>emar 
    <ranlib>emranlib 
;
EOF

    ./b2 \
        --user-config=user-config-wasm.jam \
        toolset=emscripten \
        variant=release \
        link=static \
        threading=single \
        runtime-link=static \
        --layout=system \
        --prefix="$AUX_PREFIX" \
        --without-context \
        --without-coroutine \
        --without-fiber \
        --without-thread \
        --without-stacktrace \
        --without-python \
        --without-test \
        --without-mpi \
        -j8 \
        install
)

# ------------------------------------------------------------------------------
# CDD
# ------------------------------------------------------------------------------
(
    if [[ ! -d "$EXTERN_DIR/cddlib" ]]; then
        echo "Cloning cddlib..."
        git clone https://github.com/cddlib/cddlib.git "$EXTERN_DIR/cddlib"
    fi

    cd "$EXTERN_DIR/cddlib"
    
    if [[ ! -f configure ]]; then
        ./bootstrap
    fi

    if [[ ! -f Makefile ]]; then
        CPPFLAGS="-I$AUX_PREFIX/include" \
        LDFLAGS="-L$AUX_PREFIX/lib" \
        CFLAGS="-O2" \
        CXXFLAGS="-O2" \
        emconfigure ./configure \
            --with-gmp="$AUX_PREFIX" \
            --disable-shared \
            --prefix="$AUX_PREFIX"
    fi
    
    emmake make -j8
    emmake make install

    cd "$AUX_PREFIX/include"
    ln -sf cddlib/*.h .
    ln -sf cddmp.h cdd_mp.h
)

# ------------------------------------------------------------------------------
# Frobby
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# GFAN
# ------------------------------------------------------------------------------
(
    if [[ ! -d "$EXTERN_DIR/gfanlib0" ]]; then
        echo "Downloading gfanlib..."
        mkdir -p "$EXTERN_DIR"
        curl -sSL https://users-math.au.dk/jensen/software/gfan/gfanlib0.6.2.tar.gz | tar xz -C "$EXTERN_DIR"
    fi
    
    mkdir -p "$AUX_PREFIX/include"
    ln -sfn "$AUX_PREFIX/include/cddlib" "$AUX_PREFIX/include/cdd"
    
    cd "$EXTERN_DIR/gfanlib"
    
    if [[ ! -f config.status ]]; then
        emconfigure ./configure \
            --build=i686-pc-linux-gnu \
            --host=wasm32-unknown-emscripten \
            --prefix="$AUX_PREFIX" \
            CPPFLAGS="-I$AUX_PREFIX/include -I$AUX_PREFIX/include/cddlib" \
            LDFLAGS="-L$AUX_PREFIX/lib" \
            CXXFLAGS="${CXXFLAGS:-} -std=c++11"
    fi
    
    sed -i.bak 's|../mkinstalldirs|mkdir -p|g' Makefile
    
    emmake make -j8
    emmake make install
)

# ------------------------------------------------------------------------------
# GSL
# ------------------------------------------------------------------------------
(
    mkdir -p "$AUX_BUILD/gsl"
    cd "$AUX_BUILD/gsl"
    
    if [[ ! -d "$EXTERN_DIR/gsl-2.8" ]]; then
        echo "Downloading GSL..."
        mkdir -p "$EXTERN_DIR"
        curl -sSL https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz | tar xz -C "$EXTERN_DIR"
        
        cd "$EXTERN_DIR/gsl-2.8"
        autoreconf -vfi
        cd -
    fi
    
    if [[ ! -f Makefile ]]; then
        emconfigure "$EXTERN_DIR/gsl-2.8/configure" \
            --build=i686-pc-linux-gnu \
            --host=wasm32-unknown-emscripten \
            --disable-shared \
            --prefix="$AUX_PREFIX" \
            CFLAGS="${CFLAGS:-}" \
            LDFLAGS="${LDFLAGS:-}"
    fi
    
    emmake make -j8
    emmake make install
)

# ------------------------------------------------------------------------------
# MathSAT
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# MPFR
# ------------------------------------------------------------------------------
(
    mkdir -p "$AUX_BUILD/mpfr"
    cd "$AUX_BUILD/mpfr"
    if [[ ! -d "$EXTERN_DIR/mpfr" ]]; then
        echo "Downloading MPFR source..."
        git clone https://gitlab.inria.fr/mpfr/mpfr.git "$EXTERN_DIR/mpfr"
        cd "$EXTERN_DIR/mpfr" && ./autogen.sh && cd -
    fi
    if [[ ! -f config.status ]]; then
        emconfigure "$EXTERN_DIR/mpfr/configure" \
            --build i686-pc-linux-gnu --host none \
            --with-gmp="$AUX_PREFIX" \
            --disable-shared \
            --prefix="$AUX_PREFIX"
    fi
    emmake make -j8
    emmake make install
)

# ------------------------------------------------------------------------------
# FLINT
# ------------------------------------------------------------------------------
(
    mkdir -p "$AUX_BUILD/flint"
    cd "$AUX_BUILD/flint"
    if [[ ! -d "$EXTERN_DIR/flint2" ]]; then
        echo "Cloning Flint..."
        git clone --depth 1 https://github.com/wbhart/flint2.git "$EXTERN_DIR/flint2"
        cd "$EXTERN_DIR/flint2" && ./bootstrap.sh && cd -
    fi
    if [[ ! -f Makefile ]]; then
        emconfigure "$EXTERN_DIR/flint2/configure" \
            --build=i686-pc-linux-gnu \
            --host=wasm32-unknown-emscripten \
            --with-gmp="$AUX_PREFIX" \
            --with-mpfr="$AUX_PREFIX" \
            --disable-shared \
            --disable-assembly \
            --prefix="$AUX_PREFIX"
    fi
    emmake make -j8
    emmake make install
)

# ------------------------------------------------------------------------------
# Normaliz 
# ------------------------------------------------------------------------------
(
    if [[ ! -d "$EXTERN_DIR/normaliz" ]]; then
        echo "Cloning Normaliz..."
        git clone https://github.com/Normaliz/Normaliz.git "$EXTERN_DIR/normaliz"
    fi

    cd "$EXTERN_DIR/normaliz"

    if [[ ! -f configure ]]; then
        echo "Bootstrapping Normaliz..."
        chmod +x bootstrap.sh
        ./bootstrap.sh
    fi

    if [[ ! -f Makefile ]]; then
        CPPFLAGS="-I$AUX_PREFIX/include" \
        LDFLAGS="-L$AUX_PREFIX/lib" \
        CFLAGS="-O2 -fexceptions" \
        CXXFLAGS="-O2 -fexceptions -std=c++14" \
        emconfigure ./configure \
            --build=i686-pc-linux-gnu \
            --host=wasm32-unknown-emscripten \
            --with-gmp="$AUX_PREFIX" \
            --with-flint="$AUX_PREFIX" \
            --disable-shared \
            --disable-openmp \
            --prefix="$AUX_PREFIX"
    fi
    
    emmake make -j8
    emmake make install
)

# ------------------------------------------------------------------------------
# NTL
# ------------------------------------------------------------------------------
(
    mkdir -p "$AUX_BUILD/ntl"
    cd "$AUX_BUILD/ntl"
    
    if [[ ! -d "$EXTERN_DIR/ntl" ]]; then
        echo "Cloning NTL..."
        git clone https://github.com/libntl/ntl.git "$EXTERN_DIR/ntl"
    fi
    cd "$EXTERN_DIR/ntl/src"
    
    if [[ ! -f makefile ]]; then
        emconfigure ./configure \
            CXX="em++" \
            CXXFLAGS="-O2 -fexceptions -s WASM=1 -s NODERAWFS=1" \
            PREFIX="$AUX_PREFIX" \
            GMP_PREFIX="$AUX_PREFIX" \
            NTL_GMP_LIP=on \
            NTL_STD_CXX14=on \
            SHARED=off \
            NATIVE=off \
            TUNE=generic \
            NTL_THREADS=off

        sed -e 's/^CC=gcc/CC=emcc -s NODERAWFS=1/' \
            -e 's/^WIZARD=on/WIZARD=off/' \
            makefile > makefile.patched
        mv makefile.patched makefile
    fi
    
    if ! emmake make -j8; then
        
        sed -e 's|^\t\./MakeDesc|\tchmod +x ./MakeDesc \&\& node ./MakeDesc|' \
            -e 's|^\t\./gen_gmp_aux|\tchmod +x ./gen_gmp_aux \&\& node ./gen_gmp_aux|' \
            -e 's|^\t\./gen_lip_gmp_aux|\tchmod +x ./gen_lip_gmp_aux \&\& node ./gen_lip_gmp_aux|' \
            makefile > makefile.patched
        mv makefile.patched makefile
        
        sed -i 's|if ./CheckFeatures|if node ./CheckFeatures|g' MakeCheckFeatures
        
        if [ -f MakeCheckThreads ]; then
            sed -i 's|./CheckThreads|node ./CheckThreads|g' MakeCheckThreads
        fi
        
        emmake make -j8
    fi

    emmake make install
    emranlib "$AUX_PREFIX/lib/libntl.a"
)

# ------------------------------------------------------------------------------
# CoCoA-5 Configuration & Build
# ------------------------------------------------------------------------------

if [[ -f configure ]]; then
    chmod +x configure

    sed -i.bak 's/( --with-libntl )/( --with-libntl=* )/g' configure
    sed -i.bak 's@NTL_LIB="/usr/local/lib/libntl.a".*@NTL_LIB=`echo "$option" | cut -f 2- -d=` ;;@g' configure
    sed -i.bak 's|GFAN_INC_DIR="$GFAN_LIB_DIR"|GFAN_INC_DIR="$(dirname "$GFAN_LIB_DIR")/include"|g' configure
    sed -i.bak 's|if "$PKGCONF" --exists gsl.*|if false|g' configure
    sed -i.bak 's|-lgsl-symlink  -lgslcblas  -llapack|-lgsl-symlink '"$AUX_PREFIX"'/lib/libgslcblas.a|g' configure
fi

if [[ -f configuration/boost-check-arch.sh ]]; then
    echo "#!/bin/sh" > configuration/boost-check-arch.sh
    echo "exit 0" >> configuration/boost-check-arch.sh
    chmod +x configuration/boost-check-arch.sh
fi

if [[ ! -f configuration/autoconf.mk ]]; then
    emconfigure ./configure \
        --with-cxx=em++ \
        --with-cxxflags="-O2 -fexceptions" \
        --with-libgmp="$AUX_PREFIX/lib/libgmp.a" \
        --with-boost-hdr-dir="$AUX_PREFIX/include" \
        --with-libcddgmp="$AUX_PREFIX/lib/libcddgmp.a" \
        --with-libgfan="$AUX_PREFIX/lib/libgfan.a" \
        --with-libgsl="$AUX_PREFIX/lib/libgsl.a" \
        --with-libnormaliz="$AUX_PREFIX/lib/libnormaliz.a" \
        --with-libntl="$AUX_PREFIX/lib/libntl.a" \
        --disable-mempool \
        --no-qt-gui \
        --no-readline

    ln -sf "$AUX_PREFIX/include/flint" configuration/ExternalLibs/include/flint
    ln -sf "$AUX_PREFIX/include/mpfr.h" configuration/ExternalLibs/include/mpfr.h
fi

if [[ -f configuration/autoconf.mk ]]; then
    sed -i.bak 's/-llapack//g' configuration/autoconf.mk
    sed -i.bak 's/-lblas//g' configuration/autoconf.mk
    sed -i.bak "s|-lgslcblas|$AUX_PREFIX/lib/libgslcblas.a|g" configuration/autoconf.mk
    sed -i.bak "s|-lnormaliz-symlink|-lnormaliz-symlink $AUX_PREFIX/lib/libflint.a $AUX_PREFIX/lib/libmpfr.a|g" configuration/autoconf.mk
fi

if [[ -f src/CoCoA-5/check-version-defines ]]; then
    chmod +x src/CoCoA-5/check-version-defines
fi

if [[ -f examples/CopyInfo ]]; then
    chmod +x examples/CopyInfo
fi

WASM_LDFLAGS="-O2 --fexceptions -s TOTAL_STACK=32mb -s INITIAL_MEMORY=2048mb -s ALLOW_MEMORY_GROWTH=1"

emmake make -j8 library
emmake make -j8 cocoa5 LDFLAGS="$WASM_LDFLAGS"
# emmake make -j8 -C examples executables LDFLAGS="$WASM_LDFLAGS" EXEEXT=".html"

cd src/CoCoA-5

em++ -O2 -fexceptions -s ASSERTIONS=0 -s TOTAL_STACK=32mb -s INITIAL_MEMORY=2048mb -s ALLOW_MEMORY_GROWTH=1 \
-o CoCoAInterpreter.js \
../../emscripten/wasm_patch.c \
AST.o Lexer.o Main.o Interpreter.o LineProviders.o Parser.o CoCoALibSupplement.o \
BuiltInFunctions.o BuiltInFunctions-CoCoALib.o BuiltInFunctionsVarArgs-CoCoALib.o \
BuiltInOneLiners-CoCoALib.o BuiltInFunctions-Frobby.o BuiltInFunctions-GFan.o \
BuiltInFunctions-GSL.o BuiltInFunctions-MathSAT.o BuiltInFunctions-Normaliz.o \
globals.o OnlineHelp.o VersionInfo.o Banner.o CompilationDate.o \
../../lib/libcocoa.a \
../../configuration/ExternalLibs/lib/libgfan-symlink.a \
../../configuration/ExternalLibs/lib/libcddgmp-symlink.a \
../../configuration/ExternalLibs/lib/libgsl-symlink.a \
../../extern/emscripten/install/lib/libgslcblas.a \
../../configuration/ExternalLibs/lib/libnormaliz-symlink.a \
../../extern/emscripten/install/lib/libflint.a \
../../extern/emscripten/install/lib/libmpfr.a \
../../configuration/ExternalLibs/lib/libntl-symlink.a \
../../configuration/ExternalLibs/lib/libgmpxx-symlink.a \
../../configuration/ExternalLibs/lib/libgmp-symlink.a \
../../configuration/ExternalLibs-CoCoA5/lib/libboost_filesystem-symlink.a \
--preload-file "$(pwd)/packages@/src/CoCoA-5/packages" \
--preload-file "$(pwd)/tests@/src/CoCoA-5/tests" \
--preload-file "$(pwd)/CoCoAManual@/src/CoCoA-5/CoCoAManual"

cd ../..
