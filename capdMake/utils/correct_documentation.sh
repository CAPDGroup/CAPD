#!/bin/bash -e
######################################################################################
#
#  It corrects doxygen documentation
#  - grouping: namespace capd { namespace X {  add content to group X
######################################################################################

if(($#<1))
then
  echo "USAGE $0 fileName [Y]"
  echo "   corrects documentation in the copy of fileName."
  echo "   with Y parameter applies changes to original file."
  exit 1
fi

#echo $1
FILENAME=$(basename $1)

if [[ $1 =~ ^.*/src/mpcapd/.*$ ]]
then
  arrIN=(${1//\/src\// })
  FILENAME=${arrIN[1]}
fi

#echo " filename = ${FILENAME} "

TMPFILE=tmp.file
cp $1 ${TMPFILE}
perl -i -0pe 's/^(\s*\n)+//g' ${TMPFILE}
#perl -i -0pe 's/\n\K(\s*\n)+$/\n/g' ${TMPFILE}
perl -i -0pe "s|(\@file\s*\S*)|\@file $FILENAME|" ${TMPFILE}
perl -i -0pe 's|///\s*\@addtogroup\s*\S*\s*\n///\s*@\{\s*\n||' ${TMPFILE}
perl -i -0pe 's|(namespace *capd *\{\s*\r?\n\s*namespace *(\S+) *\{ *\r?\n)|\1/// \@addtogroup \2 \n/// @\{\n|' ${TMPFILE}
perl -i -0pe 's|(namespace *capd *\{\s*\r?\n *///[^\n]*\n\s*namespace *(\S+) *\{ *\r?\n)|\1/// \@addtogroup \2 \n/// @\{\n|' ${TMPFILE}

perl -i -0pe 's|/// *\@\}\s*\n||' ${TMPFILE}
perl -i -0pe 's|(\}s*\} *// *namespace capd)|/// \@\}\n\1|' ${TMPFILE}
perl -i -0pe 's|(\}\s*\} *// *end *of *namespace capd)|/// \@\}\n\1|' ${TMPFILE}

cat ${TMPFILE}

if (( $# > 1 ))
then
  if(($2 == 'Y' )); then
    echo "APPLYING CHANGES!"
    cp ${TMPFILE} $1
  fi
fi
