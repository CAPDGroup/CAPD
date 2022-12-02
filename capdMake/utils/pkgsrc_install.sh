#!/bin/bash

set -e



install_dir=${1:-${PWD}}

echo ${install_dir}

if [ ! -e ${install_dir} ]; then
    echo "${install_dir} does not exist."
    exit 1
fi

if [[ "${install_dir}" != /* ]]
then
    echo "First argument must me an absolute path to a directory, not ${install_dir}."
    exit 1
fi


export XARGS=`which xargs`
export SH=/bin/bash



prefix_path=${install_dir}/capd_pkgsrc
export LOCAL=${prefix_path}/additional
local_src=${LOCAL}/src

pkgsrc_src=${prefix_path}/src/pkgsrc
export PKGSRC_HOME=${prefix_path}/latest
env_file=${prefix_path}/pkgsrc.env



export C_INCLUDE_PATH=${LOCAL}/include
export CPLUS_INCLUDE_PATH=${LOCAL}/include
export LIBRARY_PATH=${LOCAL}/lib



if [ ! -e $XARGS ]; then
    echo "Correct a path to xargs, (default $XARGS)"
    exit 1
fi


# pkgsrc do not like locales
unset LANG LC_ALL LANGUAGE LINGUAS
export LC_ALL=C


mkdir -p ${prefix_path}
mkdir -p ${LOCAL}
mkdir -p ${local_src}

cd ${local_src}

installed_cookie=".installed.done"

if [ ! -e termcap-1.3.1/${installed_cookie} ]; then
    rm -fr termcap-1.3.1
    wget http://ftp.gnu.org/gnu/termcap/termcap-1.3.1.tar.gz -O termcap-1.3.1.tar.gz
    tar xzf termcap-1.3.1.tar.gz
    cd termcap-1.3.1/
    ./configure --prefix ${LOCAL}
    make -j
    make  install
    touch ${installed_cookie}
fi


cd ${local_src}

if [ ! -e ncurses-5.9/${installed_cookie} ]; then
    rm -fr ncurses-5.9
    wget http://ftp.gnu.org/pub/gnu/ncurses/ncurses-5.9.tar.gz -O ncurses-5.9.tar.gz
    tar xzf ncurses-5.9.tar.gz
    cd ncurses-5.9
    mkdir -p ~/tmp
    ./configure --prefix ${LOCAL} --with-shared
    make -j
    make  install
    touch ${installed_cookie}
fi



cd ${local_src}

if [ ! -e zlib-1.2.8/${installed_cookie} ]; then
    rm -fr zlib-1.2.8
    wget "http://downloads.sourceforge.net/project/libpng/zlib/1.2.8/zlib-1.2.8.tar.gz?r=http%3A%2F%2Fzlib.net%2F&ts=1368789694&use_mirror=garr" -O zlib-1.2.8.tar.gz
    tar xzf zlib-1.2.8.tar.gz
    cd zlib-1.2.8/
    ./configure --prefix=${LOCAL}
    make -j
    make -j install
    touch ${installed_cookie}
fi


if [ ! -e ${PKGSRC_HOME}/${installed_cookie} ]; then

    if [ ! -e ${pkgsrc_src}/${installed_cookie} ]; then
	rm -fr ${pkgsrc_src}
	mkdir -p ${pkgsrc_src}/src
	cd ${prefix_path}/src
   PKGSRC_URL="http://ftp.netbsd.org/pub/pkgsrc/pkgsrc-2014Q1/pkgsrc.tar.gz"
#"ftp://ftp.netbsd.org/pub/pkgsrc/current/pkgsrc.tar.gz"
	wget ${PKGSRC_URL} -O pkgsrc.tar.gz
	tar xzf pkgsrc.tar.gz

	cd ${pkgsrc_src}

	sed -i".bac" 's|#MAKE_JOBS|MAKE_JOBS|' mk/defaults/mk.conf
	sed -i".bac" 's|/usr/include|${LOCAL}/include|' mk/term{info,cap}.builtin.mk
	sed -i".bac" 's|/usr/include|${LOCAL}/include|' mk/curses.builtin.mk
	sed -i".bac" 's|${ID} -n -G|${ID} -n -g|' mk/unprivileged.mk
	sed -i".bac" 's|#MASTER_SORT=.*|MASTER_SORT= .pl .org .de .se .no|' mk/defaults/mk.conf
	sed -i".bac" 's|#ACCEPTABLE_LICENSES=.*|ACCEPTABLE_LICENSES+=openssl|' mk/defaults/mk.conf

	touch ${pkgsrc_src}/${installed_cookie}
    fi

    cd ${pkgsrc_src}

    cd bootstrap
    mkdir -p ${PKGSRC_HOME}
    rm -fr work
    ./bootstrap --prefix ${PKGSRC_HOME} --unprivileged
    touch ${PKGSRC_HOME}/${installed_cookie}
fi

export PATH=${PKGSRC_HOME}/sbin:${PKGSRC_HOME}/bin:${PATH}
export MANPATH=${PKGSRC_HOME}/man:${MANPATH}

export C_INCLUDE_PATH=${PKGSRC_HOME}/include:${LOCAL}/include
export CPLUS_INCLUDE_PATH=${PKGSRC_HOME}/include:${LOCAL}/include
export LIBRARY_PATH=${PKGSRC_HOME}/lib:${LOCAL}/lib


bmake="${PKGSRC_HOME}/bin/bmake"

function execute_targets() {
    target="$1"
    ${bmake} ${target} -C ${pkgsrc_src}/pkgtools/pkg_install

    ${bmake} ${target} -C ${pkgsrc_src}/devel/zlib

    ${bmake} ${target} -C ${pkgsrc_src}/security/openssl

    ${bmake} ${target} -C ${pkgsrc_src}/devel/libtool-base

    ${bmake} ${target} -C ${pkgsrc_src}/devel/autoconf

    ${bmake} ${target} -C ${pkgsrc_src}/devel/automake

    ${bmake} ${target} -C ${pkgsrc_src}/devel/pkg-config

    ${bmake} ${target} -C ${pkgsrc_src}/devel/boost-libs

    ${bmake} ${target} -C ${pkgsrc_src}/math/mpfr

    ${bmake} ${target} -C ${pkgsrc_src}/x11/libX11

# ${bmake} ${target} -C ${pkgsrc_src}/graphics/tiff ALLOW_VULNERABLE_PACKAGES=yes

# ${bmake} ${target} -C ${pkgsrc_src}/x11/wxGTK28
}

echo 'PLIST_SRC+=${PKGDIR}/PLIST.common' >> ${pkgsrc_src}/security/openssl/Makefile

execute_targets install


if [ ! -e ${PKGSRC_HOME}/bin/python-config ]; then
	ln -s ${PKGSRC_HOME}/bin/python2.7-config ${PKGSRC_HOME}/bin/python-config
fi



cat > ${env_file}  <<EOF

export PATH=${PKGSRC_HOME}/sbin:${PKGSRC_HOME}/bin:\${PATH}
export MANPATH=${PKGSRC_HOME}/man:\${MANPATH}

export C_INCLUDE_PATH=${PKGSRC_HOME}/include:${LOCAL}/include:\${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${PKGSRC_HOME}/include:${LOCAL}/include:\${CPLUS_INCLUDE_PATH}
export LIBRARY_PATH=${PKGSRC_HOME}/lib:${LOCAL}/lib:\${LIBRARY_PATH}

EOF


execute_targets clean clean-depends

echo -e "!!!!!!!!!!!!!!!!!\n!!!!!!!!!!!!!!!!!\nTo use installed PKGSRC tools execute:\nsource ${env_file}"
