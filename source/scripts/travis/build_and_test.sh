#!/bin/bash
env
function or_die () {
    "$@"
    local status=$?
    if [[ $status != 0 ]] ; then
        echo ERROR $status command: $@
        exit $status
    fi
}

or_die mkdir travis-build
cd travis-build
if [[ "$DO_BUILD" == "yes" ]] ; then
    export F77=gfortran
    export CXX=`which $COMPILER`
    or_die ../configure --with-extra-cxx-flags=-std=c++11
    or_die make -j 3 VERBOSE=1
fi

exit 0
