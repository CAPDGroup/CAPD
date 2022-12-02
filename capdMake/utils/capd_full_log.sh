#!/bin/bash

CAPD_SRC_ROOT=$1

if [ -z "$CAPD_SRC_ROOT" ]; then
    echo "usage: <CAPD_SRC_ROOT>"
    exit 1
fi

CAPD_UTILS="$CAPD_SRC_ROOT/capdMake/utils"

REPOS="${CAPD_SRC_ROOT}/. "

if [ -e "${CAPD_SRC_ROOT}/capdDynSys4" ]; then
    REPOS="${REPOS} ${CAPD_SRC_ROOT}/capdDynSys4"
elif [ -e "${CAPD_SRC_ROOT}/capdDynSys" ]; then
    REPOS="${REPOS} ${CAPD_SRC_ROOT}/capdDynSys"
fi

if [ -e "${CAPD_SRC_ROOT}/capdRedHom" ]; then
    REPOS="${REPOS} ${CAPD_SRC_ROOT}/capdRedHom"
fi

if [ -e "${CAPD_SRC_ROOT}/capdExtHom" ]; then
    REPOS="${REPOS} ${CAPD_SRC_ROOT}/capdExtHom"
fi


if [ -e "${CAPD_SRC_ROOT}/.git" ]; then
    function vcs_log() {
        abs_path=`cd $1 && pwd`
        git --git-dir=$abs_path/.git --work-tree=$abs_path svn log
    }
else
    function vcs_log() {
        svn log $1
    }
fi

for repo in $REPOS; do (vcs_log $repo | ${CAPD_UTILS}/format_svn_log); done
