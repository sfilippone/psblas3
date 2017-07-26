#!/bin/sh
# $Id$

# This file is still not complete.

# this should create configure from configure.ac..

#touch NEWS AUTHORS ChangeLog COPYING INSTALL missing install-sh

#libtoolize

aclocal -I config  || { echo "no aclocal  ?" ; exit 1 ; }

# we produce a configure script
autoconf || { echo "no autoconf ?" ; exit 1 ; }

# we produce a brand new Makefile
#automake -i --add-missing || { echo "no automake ?" ; exit 1 ; }
#automake  -i || { echo "no automake ?" ; exit 1 ; }

# The automake required for autogen.sh'in this package is 1.10#.
# So users who want to build from the svn repository are required to use this version.
# 
# Users building from the tarball shouldn't bother, of course, 
# because they get the configure script generated from the tarball maintainer.

