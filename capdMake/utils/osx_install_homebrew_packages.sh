for pkg in "coreutils" "autoconf" "automake" "libtool" "pkg-config" "svn" "mpfr" "boost" "boost-python" "log4cxx" "homebrew/x11/pari" "python"; do
    brew install $pkg
done

pip install --user nose
