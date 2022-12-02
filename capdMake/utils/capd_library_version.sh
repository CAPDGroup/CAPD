#!/bin/bash

#
#
# TODO to limit number of processed commits we can start from a previous saved revision.
#
#

CAPD_SRC_ROOT="$1"
KEEP="$2"


if [ -z "${CAPD_SRC_ROOT}" ]; then
    echo "usage: <CAPD_SRC_ROOT>"
    exit 1
fi


CAPD_MAKE_PATH="${CAPD_SRC_ROOT}/capdMake"
CAPD_UTILS="${CAPD_MAKE_PATH}/utils"

VERSION_FILE="${CAPD_MAKE_PATH}/capd_version_number.m4"
LIBRARY_VERSION_FILE="${CAPD_MAKE_PATH}/capd_library_version_number.m4"
VERSION_FILE_RAW="${CAPD_MAKE_PATH}/capd_version_number.raw"

if [ -n "${KEEP}" ] && [ -e "${VERSION_FILE}" ] && [ -e "${LIBRARY_VERSION_FILE}" ] && [ -e "${VERSION_FILE_RAW}" ]; then
    exit 0
fi

echo "Regenerating version files ${VERSION_FILE} ${LIBRARY_VERSION_FILE}"

# see https://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html

current=1
revision=0
age=0

major=4
minor=1
patch=0

interface_added=no
interface_changed=no
interface_removed=no
code_changed=no

# function random_test() {
#     [ $((RANDOM % 100)) = 0 ] && commit="$commit CAPD_RELEASE:MAJOR"
#     [ $((RANDOM % 5)) = 0 ] && commit="$commit CAPD_RELEASE:MINOR"
#     [ $((RANDOM % 50)) = 0 ] && commit="$commit CAPD_INTERFACE:ADDED"
#     [ $((RANDOM % 200)) = 0 ] && commit="$commit CAPD_INTERFACE:CHANGED"
#     [ $((RANDOM % 200)) = 0 ] && commit="$commit CAPD_INTERFACE:REMOVED"
# }


function update_version() {
    if [ "${code_changed}" = "yes" ]; then
        revision=$((revision + 1))
    fi

    if [ "${interface_added}" = "yes" ] || [ "${interface_changed}" = "yes" ] || [ "${interface_removed}" = "yes" ]; then
        revision=0
        current=$((current + 1))
    fi

    if [ "${interface_added}" = "yes" ]; then
        age=$((age + 1))
    fi

    if [ "${interface_changed}" = "yes" ] || [ "${interface_removed}" = "yes" ]; then
        age=0
    fi
}

function update_changelog() {
    changelog_file=${CAPD_MAKE_PATH}/release_changelog/ChangeLog-${major}.${minor}
    if [ ! -e ${changelog_file} ]; then
        echo "ERROR: Cannot find release ChangeLog file: ${changelog_file}" >& 2
        
				
    else

        changelog_header=$(head -n 1 ${changelog_file})
        if ! echo "${changelog_header}" | grep "capd (${major}.${minor}).*" ; then
            echo "ERROR: Cannot find correct release version (capd (${major}.${minor})) in ChangeLog file header: ${changelog_file} ${changelog_header}" >& 2
          
        else 
            cat ${changelog_file} | cat - ${CAPD_SRC_ROOT}/ChangeLog > ${CAPD_SRC_ROOT}/ChangeLog.tmp && mv ${CAPD_SRC_ROOT}/ChangeLog.tmp ${CAPD_SRC_ROOT}/ChangeLog
        fi
		fi
}

echo > ${CAPD_SRC_ROOT}/ChangeLog
prev_date=""

while read log ; do
    commit=$(echo $log | cut -d'|' -f4)
    date_time=$(echo $log | cut -d'|' -f1)
    date=$(echo ${date_time} | sed 's/\(.*\)_.*/\1/')
#    echo $date_time $commit
#    random_test

    release_changed=no
    code_changed=yes

    if [ "$date" != "$prev_date" ]; then
        patch=$((patch + 1))
    fi

    if echo ${commit} | grep -q "CAPD_RELEASE:MAJOR"; then
        major=$((major + 1))
        minor=0
        patch=0
        release_changed=yes
    fi

    if echo ${commit} | grep -q "CAPD_RELEASE:MINOR"; then
        minor=$((minor + 1))
        patch=0
        release_changed=yes
    fi


    if echo ${commit} | grep -q "CAPD_INTERFACE:ADDED"; then
        interface_added=yes
    fi

    if echo ${commit} | grep -q "CAPD_INTERFACE:CHANGED"; then
        interface_changed=yes
    fi

    if echo ${commit} | grep -q "CAPD_INTERFACE:REMOVED"; then
        interface_removed=yes
    fi

    if [ "${release_changed}" = "yes" ]; then
        update_version
#        update_changelog

        interface_added=no
        interface_changed=no
        interface_removed=no
        code_changed=no

        echo "${major}.${minor}.${patch} ${current}:${revision}:${age}"
    fi

    prev_date="${date}"
done < <(${CAPD_UTILS}/capd_full_log.sh ${CAPD_SRC_ROOT} | sort)

update_version

#it is required by debian package builder to have exact last version number
#echo -en "capd (${major}.${minor}.${patch}-1) stable; urgency=low\n\n  * Automatic ChangeLog. Patch ${major}.${minor}.${patch}\n\n -- Automatic <no_email@noemail.com>  $(date +'%a, %d %b %Y %T %z')\n\n" | cat - ${CAPD_SRC_ROOT}/ChangeLog > ${CAPD_SRC_ROOT}/ChangeLog.tmp && mv ${CAPD_SRC_ROOT}/ChangeLog.tmp ${CAPD_SRC_ROOT}/ChangeLog


echo "Saving ${major}.${minor}.${patch} ${current}:${revision}:${age}"

echo -e "#file genereted using capd_library_version.sh\nm4_define([CAPD_VERSION_NUMBER], [${major}.${minor}.${patch}])" > "${VERSION_FILE}"

echo -e "#file genereted using capd_library_version.sh\nm4_define([CAPD_LIBRARY_VERSION_NUMBER], [${current}:${revision}:${age}])" >  "${LIBRARY_VERSION_FILE}"

echo -e "${major}.${minor}.${patch}" > "${VERSION_FILE_RAW}"
