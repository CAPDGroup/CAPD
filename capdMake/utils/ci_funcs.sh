#!/bin/bash


function go_to_dist() {
    files=($PWD/input/capd-*[0-9]*.tar.gz)

    if [ "1" != "${#files[*]}" ]; then
        echo "ERROR: Found more than one file"
        exit 1
    fi

    mkdir workdir
    cd workdir

    dist_archive=${files[0]}
    tar xzf ${dist_archive}

    dirs=(capd-*/)

    if [ "1" != "${#dirs[*]}" ]; then
        echo "ERROR: Found more than one directory"
        exit 1
    fi

    dist=${dirs[0]}
    cd $dist
}

function pwd_version() {
    echo $(basename $PWD) | sed 's/capd-\(.*\)/\1/'
}

function detect_cores() {
    if [ -z "$CORES" ]; then
        if [ $(uname) == 'Darwin' ]; then
            CORES=$(sysctl hw.ncpu | awk '{print $2}')
        elif [ -e /proc/cpuinfo ]; then
            CORES=$(grep -c ^processor /proc/cpuinfo)
        else
            CORES="1"
        fi
        echo "Using $CORES CPU cores"
    else
        echo "CORES set to $CORES"
    fi
}

function do_distcheck() {
    ppwd=$(pwd)

    BUILD_DIR="$1"
    MAKE_ARGS="$2"

    detect_cores

#    autoreconf --force --install --verbose


    if [ -e $BUILD_DIR ]; then
        chmod -R u+rw $BUILD_DIR
        rm -r $BUILD_DIR
    fi

    mkdir -p $BUILD_DIR
    cd $BUILD_DIR

   if [ -n "${CROSS_TARGET_HOST}" ]; then
      echo "Cross platform configure ${CROSS_TARGET_HOST}"
      ../configure --host=${CROSS_TARGET_HOST} --prefix $PWD/install/
   else
      echo "Default configure"
      ../configure --prefix $PWD/install/
   fi

   make -j $CORES distcheck $MAKE_ARGS

   cd $ppwd
}


function do_dist() {
    ppwd=$(pwd)

    BUILD_DIR="$1"
    MAKE_ARGS="$2"
    OUTPUT_DIR="$3"

    if [ -e $BUILD_DIR ]; then
        chmod -R u+rw $BUILD_DIR
        rm -r $BUILD_DIR
    fi

    mkdir -p $BUILD_DIR
    cd $BUILD_DIR

    ../configure

    if command -v zip 2>/dev/null; then
        make -j $CORES dist DIST_TARGETS="dist-gzip dist-zip" $MAKE_ARGS
    else
        echo "Do not have zip"
        make -j $CORES dist DIST_TARGETS=dist-gzip $MAKE_ARGS
    fi

    if [ ! -z "$OUTPUT_DIR" ]; then
        mkdir -p "$OUTPUT_DIR"
        mv *.{tar.gz,zip} $OUTPUT_DIR/
    fi;

    cd $ppwd
}
