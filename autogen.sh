#!/bin/sh

aclocal -I m4
if which libtoolize >/dev/null; then
  libtoolize --copy
else
  if which glibtoolize >/dev/null; then
    glibtoolize --copy
  fi
fi
#libtoolize --copy
autoconf
autoheader
automake --add-missing --copy
