#!/bin/sh

dir=$1;

dir_cmp=`echo $dir | sed 's./. /.g'`
if [ ! -d $dir ]
then
  path=''
  for cmp in $dir_cmp ;   do
      path="$path$cmp";
      if [ ! -d $path ] ; then 
	 mkdir $path; rc=$?;
	 if [ $rc != 0 ] ; then 
	    echo "Error while making directory $path "
	    exit 1
	 fi
      fi
  done
fi
        