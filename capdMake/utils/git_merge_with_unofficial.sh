#!/bin/bash

FROM_BRANCH='master_unofficial'


for file in $(git diff --name-only ${FROM_BRANCH} | grep -v multiVectorField ); do
    if git cat-file -e ${FROM_BRANCH}:${file} &> /dev/null ; then
        echo "File ${file} exists on the branch"
        git checkout ${FROM_BRANCH} $file
    else
        git rm $file
        echo "File ${file} does not exists on the branch"
    fi

done
