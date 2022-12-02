version="2.7.11"

curl  -L -o /tmp/XQuartz-${version}.dmg "https://dl.bintray.com/xquartz/downloads/XQuartz-${version}.dmg"
hdiutil attach /tmp/XQuartz-${version}.dmg
sudo installer -verbose -pkg /Volumes/XQuartz-${version}/XQuartz.pkg -target /
ln -fs /opt/X11/include/X11 /usr/local/include/X11
ln -fs /opt/X11/lib/libX11.dylib /usr/local/lib/libX11.dylib
ln -fs /opt/X11/lib/pkgconfig/{x11,xproto,kbproto,xcb,pthread-stubs,xau,xdmcp}.pc /usr/local/lib/pkgconfig/
