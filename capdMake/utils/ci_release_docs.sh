#!/bin/bash

set -e

export PATH=/usr/local/bin:$PATH
source "$(dirname $0)/ci_configure_flags.sh"
source "$(dirname $0)/ci_funcs.sh"

env


go_to_dist

if ! [ -e capdDynSys4 ]; then
  ln -s capdDynSys capdDynSys4 # it is required in Doxyfile
fi

if ! [ -e capdDynSys ]; then
  ln -s capdDynSys4 capdDynSys
fi

(cd capdMake/docs && doxygen Doxyfile)

(cd capdDynSys/docs && make doc)
tar czf capdDynSys-docs.tar.gz capdDynSys/docs/html

(cd capdRedHom/docs && make doc)
tar czf capdRedHom-docs.tar.gz capdRedHom/docs/html
