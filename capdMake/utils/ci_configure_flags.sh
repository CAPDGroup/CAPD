#!/bin/sh

# set -o verbose

#
# DO WE NEDD THIS?? MJ
#

# # it is required by nested configure in distcheck because we cannot pass there any arguments
# export with_boost=${WITH_BOOST}
# boost_version_req="104900"
# tmpfile="$(mktemp XXXXXXX).cpp"
# cat - <<_ACEOF > ${tmpfile}
# #include <boost/version.hpp>

# #if !defined BOOST_VERSION
# # error BOOST_VERSION is not defined
# #elif BOOST_VERSION < $boost_version_req
# # error Boost headers version < $boost_version_req
# #endif

# int main ()
# {
#   return 0;
# }
# _ACEOF

# # export with_boost, so we can use it with source
# if g++ -c "${tmpfile}" -o "${tmpfile}.o" > /dev/null 2>&1 ; then
#     export with_boost="yes"
#     echo "Boost detected"
# else
#     export with_boost="no"
#     echo "Boost not detected"
# fi

# rm -f "${tmpfile}" "${tmpfile}.o"
