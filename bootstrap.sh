#!/bin/bash
autoheader
aclocal -I m4 --install
libtoolize --automake
automake --add-missing
autoconf
