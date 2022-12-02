#!/bin/bash

set -e


#it is quite difficult to set PATH for jenkins-slave started by launchd on OSX.
export PATH=/usr/local/bin:$PATH
source "$(dirname $0)/ci_configure_flags.sh"
source "$(dirname $0)/ci_funcs.sh"


output=$PWD/output
input=$PWD/input

go_to_dist

./bootstrap.sh

ver="unofficial-$(pwd_version)"
do_dist build "VERSION=${ver}" "$output"
mv $input/build_date $output/
