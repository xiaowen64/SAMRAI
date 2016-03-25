#!/bin/sh

##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/testtools/testcount.sh $
## Package:     SAMRAI test
## Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedVersion$
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: script for testing the number of tests passed
##

num_tests=`echo $1 |  sed -e 's/,/ /g' | wc -w`
num_tests=`expr $num_tests \* $2`
num_passed=`grep PASSED $3 | wc -l`
if test $num_tests -ne $num_passed; 
   then echo "FAILED:  Count of passed tests in directory " `pwd` " passed was $num_passed expected $num_tests"
fi
