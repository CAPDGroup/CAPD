#!/bin/bash

set -e

if ! [ -z "${TARGET_HOST}" ]; then
    echo "TARGET_HOST=${TARGET_HOST}"
else

    if [ -z "$NODE_NAME" ]; then
        NODE_NAME="jenkins-slave-noname"
    fi

    case "${NODE_NAME}" in
        jenkins-slave-noname)
            TARGET_HOST="noname" ;;
        *-jenkins-slave-dockerin-*)
            TARGET_HOST="$(echo ${docker_image} | sed 's/.*-\(.*\)/\1/')" ;;
        *jenkins-slave-*)
            TARGET_HOST="$(echo $NODE_NAME  | sed 's/.*jenkins-slave-\([a-zA-Z0-9._]*\)-[0-9a-zA-Z]*$/\1/')";;

    esac

    if [ -z "${TARGET_HOST}" ]; then
        echo "Cannot guess TARGET_HOST from ${NODE_NAME}" >& 2
        exit 1
    fi

    echo "TARGET_HOST=${TARGET_HOST}"
    export TARGET_HOST
fi
