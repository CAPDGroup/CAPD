#!/bin/bash

set -e

#it is quite difficult to set PATH for jenkins-slave started by launchd on OSX.
export PATH=/usr/local/bin:$PATH
source "$(dirname $0)/ci_configure_flags.sh"
source "$(dirname $0)/ci_funcs.sh"
source "$(dirname $0)/ci_target_host.sh"
export WITHOUT_CAPD_EXAMPLES=true

env

root_dir=$PWD

go_to_dist

mkdir -p build_result
./configure

make -C capdRedHom/programs/apiRedHom-py py_src
cp capdRedHom/programs/apiRedHom-py/*.tar.gz build_result
