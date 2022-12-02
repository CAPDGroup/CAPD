#!/bin/bash

set -e

if [ -e capdRedHom ]; then
    (cd capdRedHom/ && ./bootstrap.sh)
fi

autoreconf --force --install --verbose
