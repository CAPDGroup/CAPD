#!/bin/bash


set -e

#it is quite difficult to set PATH for jenkins-slave started by launchd on OSX.
export PATH=/usr/local/bin:$PATH
source "$(dirname $0)/ci_configure_flags.sh"
export WITHOUT_CAPD_EXAMPLES=true

env


files=(*.tar.gz)

if [ "1" != "${#files[*]}" ]; then
   echo "ERROR: Found more than one file"
   exit 1
fi

dist_archive=${files[0]}
version=$(echo $dist_archive | sed 's/capd.*-\(.*\)\.tar\.gz/\1/')

cp ./capdMake/osx/capd.rb.in capd.rb
sed -i.bac "s/src_SHA1_PLACEHOLDER/$(gsha1sum ${dist_archive}  | cut -f1 -d' ')/" capd.rb
sed -i.bac "s|src_URL_PLACEHOLDER|file://${PWD}/${dist_archive}|" capd.rb
sed -i.bac "s|VERSION_PLACEHOLDER|${version}|" capd.rb

rm -f /Library/Caches/Homebrew/capd* /Library/Caches/Homebrew/Formula/capd*
brew uninstall --force capd
brew install --build-bottle capd.rb
brew bottle capd.rb
brew uninstall --force capd
rm -f /Library/Caches/Homebrew/capd* /Library/Caches/Homebrew/Formula/capd*
