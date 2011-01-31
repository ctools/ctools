#!/bin/sh

aclocal -I m4
libtoolize --copy
autoconf
autoheader
automake --add-missing --copy
