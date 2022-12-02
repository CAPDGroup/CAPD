#!/bin/bash
######################################################################################
#
#  It generates Makefile.am containing list of all header files in current directory
#
######################################################################################

noinst=$1

if [ "$(uname)" == 'Darwin' ]; then
    SED='gsed'
else
    SED='sed'
fi


echo 'include ${capdMake}/make/common_makefile.mkf'
echo

if [ -z $noinst ]; then
    echo 'nobase_include_HEADERS = \'
else
    echo 'noinst_HEADERS = \'
fi

find . \( -name "*.h" -o -name "*.hpp" \)| sort | ${SED} -n '$!{s/\(.*\)/\1 \\/};s/^\.\/\(.*\)/\1/;p;'
echo
echo
