#!/bin/bash

set -e

function do_config() {
   if [ -n "${CROSS_TARGET_HOST}" ]; then
      echo "Cross platform ./configure --host=${CROSS_TARGET_HOST}"
      ./configure --host=${CROSS_TARGET_HOST}
   else
      echo "Default configure"
      ./configure
   fi
}



function do_relese_build() {
    ppwd=$(pwd)
    output="$1"

    ver="$(pwd_version)"

    mkdir -p "$output"

    if [ $(uname) == 'Darwin' ]; then
        CORES=$(sysctl hw.ncpu | awk '{print $2}')
    elif [ -e /proc/cpuinfo ]; then
        CORES=$(grep -c ^processor /proc/cpuinfo)
    else
        CORES=""
    fi

    echo "Detected $CORES CPU cores"

    if [ -z "${TARGET_HOST}" ]; then
        echo "Please set TARGET_HOST, e.g. ubuntu_14.04" >& 2
        exit 1
    fi

    echo "TARGET_HOST=${TARGET_HOST}"

    mkdir build_result
    do_config

    if [ -n "${CROSS_TARGET_HOST}" ] && [ "${CROSS_TARGET_HOST}" != "i686-pc-linux-gnu" ] ; then
        echo "Only build for cross-compilation"
        # run a few times, to avoid infrastructure errors during the build
        (make -j $CORES) || (make -j $CORES) || (make -j $CORES)
    else
        echo "Build and check"
        # run a few times, to avoid infrastructure errors during the build        
        (make -j $CORES check) || (make -j $CORES check) || (make -j $CORES check)
    fi
    
    if [ -n "${CROSS_TARGET_HOST}" ]; then
        BIN_NAME="capd-${ver}-dev-${CROSS_TARGET_HOST}"
    else
        BIN_NAME="capd-${ver}-dev-${TARGET_HOST}-$(uname -m)"
    fi

    make install DESTDIR=$PWD/$BIN_NAME

   zip -r $BIN_NAME.zip $BIN_NAME
   mv $BIN_NAME.zip $output/
   tar zcvf $BIN_NAME.tar.gz $BIN_NAME
   mv $BIN_NAME.tar.gz $output/

   cd $ppwd
}

#it is quite difficult to set PATH for jenkins-slave started by launchd on OSX.
# export PATH=/usr/local/bin:$PATH

source "$(dirname $0)/ci_configure_flags.sh"
source "$(dirname $0)/ci_funcs.sh"
source "$(dirname $0)/ci_target_host.sh"

export WITHOUT_CAPD_EXAMPLES=true


output=$PWD/output
input=$PWD/input

env


go_to_dist
do_relese_build "$output"
mv $input/build_date $output/
