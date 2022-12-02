#!/bin/bash


makefiles="$(find -L . -name Makefile.am)"

for file in $(find -L . -name "*.h" -o -name "*.hpp" -o -name "*.c" -o -name "*.cpp" -o -name "*.py" -o -name "*.png" -o -name "*.jpg"); do

    file_bn=$(basename $file)
    if ! grep --word-regexp -q "$file_bn" $makefiles; then
        echo "Cannot find $file in any Makefile.am"
    fi
done
