#!/bin/bash

set -e


src_file=(src/capd*.tar.gz)
src_sha=$(sha1sum ${src_file}  | cut -f1 -d' ')
version=$(echo $src_file | sed 's/.*capd.*-\(.*\)\.tar\.gz/\1/')

sf_url="http://sourceforge.net/projects/capd/files/${version}"

mkdir -p python bin src

cat - > bin/README.md <<EOF
This folder contains compiled CAPD.

*  *amd64.deb - Ubuntu package compiled on Ubuntu 13.10
*  capd-*-dev-<distro>.tar.gz - binaries compiled on linux <distro> with prefix=/usr/local
*  capd-capdDynSys*-dev-<distro>.tar.gz - binaries of capdDynSys compiled on linux <distro> with prefix=/usr/local
*  capd-*-dev-osx-<version>.tar.gz - binaries compiled on Mac OSX <version>
*  capd-*-<version>.bottle.tar.gz - binaries compiled on Mac OSX <version> for Homebrew
*  *mingw*.zip - cross-compiled binaries for Windows using Linux [MXE](http://mxe.cc/)
*  capd.rb - formula for OSX [Homebrew](http://brew.sh/) package manager, to instale type: brew install ${sf_url}/bin/capd.rb
*  capdRedHomM-<version>-<distro>.zip - Mathematica bindings for <distro>.

Let us know if you need more,
Enjoy!

EOF

cat - > src/README.md <<EOF
This folder contains CAPD source files.
* capd-*.tar.gz - sources of capdDynSys and capdRedHom
* capd-capdDynSys*.tar.gz - sources of capdDynSys
Enjoy!

EOF

cat - > README.md <<EOF
This folder contains CAPD files.
 * src - source code.
 * bin - compiled.
Enjoy!

EOF


#mv -n capd-release-osx/*/* bin/
#mv -n capd-release-osx-bottle/*/* bin/
mv -n capd-release-bin/* bin/

# for deb_dir in capd-release-deb-docker/*; do
#     system=$(echo $deb_dir | sed 's/.*docker_image=capd-api-\(.*\),.*/\1/')
#     mkdir -p bin/${system}
#     mv -n ${deb_dir}/* bin/${system}/
# done


# cp ./capdMake/osx/capd.rb.in bin/capd.rb

# echo "Replace homebrew src variables"
# sed -i.bac "s|src_SHA1_PLACEHOLDER|${src_sha}|" bin/capd.rb
# sed -i.bac "s|src_URL_PLACEHOLDER|${sf_url}/${src_file}|" bin/capd.rb
# sed -i.bac "s|VERSION_PLACEHOLDER|${version}|" bin/capd.rb
# sed -i.bac "s|bottle_URL_PLACEHOLDER|${sf_url}/bin|" bin/capd.rb


# for bottle_file in bin/capd*-${version}.*.bottle.tar.gz; do
#     echo "Replace homebrew bottle variables ${bottle_file}"
#     osx_version_name=$(basename ${bottle_file} | sed "s|capd.*-${version}\.\(.*\)\.bottle\.tar\.gz|\1|")
#     echo ${osx_version_name}
#     bottle_sha=$(sha1sum ${bottle_file}  | cut -f1 -d' ')
#     sed -i.bac "s|${osx_version_name}_SHA1_PLACEHOLDER|${bottle_sha}|" bin/capd.rb
# done

# rm bin/capd.rb.bac




#setup_ssh_keys
cat - > ssh_askpass_env <<EOF
#!/bin/sh
echo "\$KEYS_PASSPHRASE"
EOF

chmod a+x ssh_askpass_env

export DISPLAY=0.0
export SSH_ASKPASS=$PWD/ssh_askpass_env
/usr/bin/keychain --timeout 60 $HOME/.ssh/id_rsa
source $HOME/.keychain/$HOSTNAME-sh

echo "SSH add"
ssh-add -l

echo "rsync"

find src bin  python -name ".*" | xargs rm -fr
rsync -aivr src bin  python README.md capd,capd@frs.sourceforge.net:/home/frs/project/capd/${version}
echo -e "rm newest \n ln -s ${version} newest" | sftp  capd,capd@frs.sourceforge.net:/home/frs/project/capd


echo "build_www"

rm -fr build_www
mkdir build_www
cd build_www
tar xzf ../www/capdWWW.tar.gz

echo "rsync www"
rsync -aivr . capd,capd@web.sourceforge.net:/home/project-web/capd/htdocs/
cd -

echo "kill agent"
kill $SSH_AGENT_PID

echo "pypi_pass"
pypi_pass=$(echo $ENCODED_PYPI_PASS | openssl enc -aes-256-cbc  -a -d -pass env:KEYS_PASSPHRASE)

echo "twin upload"
# twine upload  -u mateusz.juda python/capdRedHomPy-${version}.tar.gz -p ${pypi_pass}
