#!/bin/bash

set -e

SUB_PACKAGE="$1"

#it is quite difficult to set PATH for jenkins-slave started by launchd on OSX.
export PATH=/usr/local/bin:$PATH
source "$(dirname $0)/ci_configure_flags.sh"
source "$(dirname $0)/ci_funcs.sh"

env

go_to_dist

if [ "$SUB_PACKAGE" == "capdDynSys" ]; then
    rm -fr capdRedHom
    rm -fr capdExtHom
elif [ "$SUB_PACKAGE" == "capdRedHom" ]; then
    rm -fr capdDynSys
else
   echo "Wrong package ${SUB_PACKAGE}" > 2
   exit 1
fi

do_distcheck "build_snapshot_${SUB_PACKAGE}" "PACKAGE=${SUB_PACKAGE} VERSION=$(pwd_version)"
