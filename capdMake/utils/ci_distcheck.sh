#!/bin/bash

set -e

#it is quite difficult to set PATH for jenkins-slave started by launchd on OSX.
export PATH=/usr/local/bin:$PATH
source "$(dirname $0)/ci_configure_flags.sh"
source "$(dirname $0)/ci_funcs.sh"

export CXXFLAGS="-s"
env


go_to_dist
do_distcheck build_distcheck ""
