#!/bin/bash
set -e

docker_image="$1"
cmd="$2"

export docker_image
source "$(dirname $0)/ci_target_host.sh"

cid_file=$(mktemp)
rm ${cid_file}


volumes=""
user_local_dir="/mnt/usr_local"

[ -e "${user_local_dir}/${TARGET_HOST}/bin" ] && volumes="$volumes -v ${user_local_dir}/${TARGET_HOST}/bin:/usr/local/bin"
[ -e "${user_local_dir}/${TARGET_HOST}/lib/MATLAB" ] && volumes="$volumes -v ${user_local_dir}/${TARGET_HOST}/lib/MATLAB:/usr/local/lib/MATLAB -e MATLAB_PREFDIR=/usr/local/lib/MATLAB/R2014b" # somehow matlab cannot guess PREFDIR from mounted volume

[ -e "${user_local_dir}/${TARGET_HOST}/lib/Wolfram" ] && volumes="$volumes  -v ${user_local_dir}/${TARGET_HOST}/lib/Wolfram:/usr/local/lib/Wolfram"


docker -H unix:///run/docker.sock pull capd/"${docker_image}" # refresh

# network required for svn. We use svn command to generate next version number of the library
docker -H unix:///run/docker.sock run --network vnet --rm --cidfile="${cid_file}" --user="$(id -u)" --workdir="$PWD" --volumes-from="$(hostname -s)" $volumes -e "TARGET_HOST=${TARGET_HOST}" -e "CAPD_MANUAL_VERSION=${CAPD_MANUAL_VERSION}"  capd/"${docker_image}" /bin/bash -c "${cmd}"

#cid=$(cat "${cid_file}")
#docker stop "${cid}"
#docker rm "${cid}"
