dnl $Id$

dnl Define macros for supporting conduit.


AC_DEFUN([CASC_SUPPORT_CONDUIT],[
dnl Support conduit libraries by setting the variables
dnl conduit_PREFIX, conduit_INCLUDES, and conduit_LIBS.
dnl Arg1: empty if you want the default to be off.
dnl
# Begin macro CASC_SUPPORT_CONDUIT
CASC_ARG_WITH_ENV_WRAPPER(conduit, conduit_PREFIX,
ifelse($1,,
[  --with-conduit[=PATH]	Use conduit and optionally specify where it is installed.],
[  --without-conduit	Do not use the conduit library.]),
ifelse($1,,if test "$with_conduit" = '' ; then unset conduit_PREFIX; else conduit_PREFIX=; fi, conduit_PREFIX=)
)
CASC_VAR_SET_CONDUIT(conduit_PREFIX,conduit_INCLUDES,conduit_LIBS)
CASC_AC_LOG_VAR(conduit_PREFIX conduit_INCLUDES conduit_LIBS)
if test "${conduit_PREFIX+set}" = set; then
  btng_save_cppflags=$CPPFLAGS

  # Add conduit include flags to cpp so we can examine its header file.
  CPPFLAGS="$conduit_INCLUDES $CPPFLAGS"
  CASC_AC_LOG_VAR(conduit_INCLUDES CPPFLAGS)

  # Check if conduit header is ok.
  AC_CHECK_HEADER(conduit.hpp,:,
    [AC_MSG_ERROR(Problems checking conduit.hpp)])

  # Reset cpp after checking conduit header file.
  CPPFLAGS=$btng_save_cppflags
  unset btng_save_cppflags

  CASC_AC_LOG_VAR(CPPFLAGS)
fi
# End macro CASC_SUPPORT_CONDUIT
])


AC_DEFUN([CASC_VAR_SET_CONDUIT],[
dnl Provides support for the conduit libraries.
dnl
dnl Arguments are:
dnl 1. Name of variable to set to path where conduit is installed.
dnl    Nothig is done if this variable is unset.
dnl 2. Name of the INCLUDES variable similar to the automake INCLUDES variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl 3. Name of the LIBS variable similar to the automake LIBS variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl
# Begin macro CASC_VAR_SET_CONDUIT
if test "${$1+set}" = set ; then
  if test ! "${$2+set}" = set ; then
    test -n "${$1}" && $2="-I${$1}/include/conduit"
  fi
  if test ! "${$3+set}" = set ; then
    $3='-lconduit'
    if test -n "${$1}" ; then
      for i in ${$3} ; do
	tmp_name=`echo $i | sed 's/^-l//'`
        if test ! -f "${$1}/lib/lib${tmp_name}.a" && \
          test ! -f "${$1}/lib/lib${tmp_name}.so"; then
          AC_MSG_WARN(Library file for ${tmp_name} is missing from ${$1}/lib.)
        fi
      done
      $3="-L${$1}/lib ${$3}"
    fi
  fi
fi
# End macro CASC_VAR_SET_CONDUIT
])dnl
