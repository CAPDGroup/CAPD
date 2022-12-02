#!/bin/bash

set -e

#it is quite difficult to set PATH for jenkins-slave started by launchd on OSX.
export PATH=/usr/local/bin:$PATH

source "$(dirname $0)/ci_configure_flags.sh"

env

    if [ $(uname) == 'Darwin' ]; then
        CORES=$(sysctl hw.ncpu | awk '{print $2}')
    elif [ -e /proc/cpuinfo ]; then
        CORES=$(grep -c ^processor /proc/cpuinfo)
    else
        CORES=""
    fi

    echo "Detected $CORES CPU cores"

function do_config() {
   rm -fr lib/* bin/*
   ./configure --prefix $PWD/install/ CXXFLAGS="-Os -s" || ( autoreconf -fi && ./configure --prefix $PWD/install/ CXXFLAGS="-Os -s" ) #optimize size because Debian Whezzy linker cannot handle big files
}

function do_make() {
   JOB="$1"

   echo "do_make ${JOB}"
   make -j ${CORES} ${JOB} || ( do_config && (make -k distclean || true) && do_config && make -j ${CORES} ${JOB} ) || ( rm -fr configure capd.cache && do_config && make -j ${CORES} ${JOB})
}


do_make check
