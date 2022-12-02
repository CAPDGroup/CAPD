#!/bin/bash -xv

########################################################################
#
# USAGE:
#     create_snapshoot  [zip_name] [packages] [root_dir] [svn_dir]
#
#  - checks out main capd repository and repositories containing packages
#      into svn_dir,
#  - creates distribution (stripping SVN info) in root_dir
#  - zips it into zip_name
#  - makes compilation tests
#  - publish it to WWW directory if tests succeeded
#
#
#  NOTE:
#    All directories should be relative path
#    root_dir     need to be direct subdirectory of current directory
#
#
########################################################################


####
# zip (snapshoot) file name
zip_name=$1
: ${zip_name:="capd_src.zip"}

# modules to be included
packages=$2
: ${packages:="capdDynSys4"}
#: ${packages:="capdDynSys capdRedHom capdExtHom"}

#  name of the root directory (in the zip file)
root_dir=$3
: ${root_dir:="capd"}

# working directory
work_dir=$4
 : ${work_dir:="${PWD}/capd_snapshoot"}
mkdir -p ${work_dir}

# directory for testing compilation system of CAPD
test_dir="capd_test_dir"
mkdir -p ${test_dir}

# www directory - where to publish results
www_dir="${PWD}/www"
mkdir -p ${www_dir}

svn_cmd="svn --non-interactive --trust-server-cert "

# svn connection (choose one of the following)
svn_url=https://svn.capdnet.ii.uj.edu.pl/
#svn_url=file:///capd/

#user=kapela
#svn_url=svn+ssh://${user}@mnich.ii.uj.edu.pl/capd/

################ end of settings #######################


cd ${work_dir}


# cheking out CAPD from repository
rm -rf ${root_dir}
${svn_cmd} co ${svn_url}capd  ${root_dir}
cd ${root_dir}
for module in ${packages}
do
${svn_cmd} co ${svn_url}${module}
done

echo "generating Makefiles "
cd capdMake/utils
./generateMakefiles.sh

echo "generating configure scripts"
cd ../../
autoreconf -fi

#echo "recursively removing .svn folders from"
#pwd
#rm -rf `find . -type d -name .svn`

echo "making zip file without SVN informations"
cd ..
rm -f ${zip_name}
find ${root_dir} -type d -name .svn | sed "s/\(.*\)/\1\/\n\1\/*/" > list.ex
zip -u -9 -r ${zip_name} ${root_dir} -x@list.ex
chmod 664 ${zip_name}

echo "making compilation tests"
rm -rf ${test_dir}
mkdir -p ${test_dir}
#cp ${zip_name} ${test_dir}/
unzip ${zip_name} -d ${test_dir}
cd ${test_dir}
../${root_dir}/capdMake/utils/check_capd.sh ${root_dir}

if [[ $? -ne 0 ]] ; then
  # email subject
  SUBJECT="Compilation error"
  # Email To ?
  EMAIL="kapela@ii.uj.edu.pl"
  # Email text/message
  EMAILMESSAGE="./emailmessage.txt"
  echo "CAPD compilation error:"> $EMAILMESSAGE
  echo " zip name : ${zip_name} \n packages ${packages} \n root dir : ${root_dir}\n" >>$EMAILMESSAGE
  # send an email using /bin/mail
  /usr/bin/mail -s "$SUBJECT" "$EMAIL" < $EMAILMESSAGE
  echo "CAPD compilation error! Aborting..."
  exit 1;
fi

cd ${work_dir}

echo "publishing to WWW"
#mv --backup=numbered ${zip_name} ${www_dir}
mv  ${zip_name} ${www_dir}
#rm -f /srv/www/htdocs/capd/${snapshoot_name}.~10~
