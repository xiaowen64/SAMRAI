dnl Define a macro for supporting SILO

AC_DEFUN([CASC_SUPPORT_SILO],[

# Begin CASC_SUPPORT_SILO
# Defines SILO_PREFIX SILO_INCLUDES and SILO_LIBS if with-silo is specified.
AC_ARG_WITH(silo,
[ --with-silo=PATH  Use SILO and optionally specify where SILO is installed.],
, with_silo=no)

BTNG_AC_LOG(with_silo is $with_silo)

case "$with_silo" in
  no)
    BTNG_AC_LOG(Not setting up for SILO)
    : Do nothing
  ;;
  yes)
    # SILO install path was not specified.
    # Look in a couple of standard locations to probe if 
    # SILO header files are there.
    BTNG_AC_LOG(Looking for SILO installation)
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/silo.h; then
        silo_PREFIX=${dir}
	BTNG_AC_LOG(Found silo.h ${dir}/include)
        break
      fi
    done
  ;;
  *)
    # SILO install path was specified.
    BTNG_AC_LOG(Expect SILO installation in $with_silo)
    silo_PREFIX=$with_silo
    silo_INCLUDES="-I${silo_PREFIX}/include"
    BTNG_AC_LOG(Set silo_INCLUDES to $silo_INCLUDES)
  ;;
esac

# Use HDF silo library if available
if test "${silo_PREFIX+set}" = set; then
   if test -f ${silo_PREFIX}/lib/libsilo.a; then
      silo_LIBS='-lsilo'
   elif test -f ${silo_PREFIX}/lib/libsiloh5.a; then
      silo_LIBS='-lsiloh5'
   else
      AC_MSG_ERROR([Could not fine silo library in $silo_PREFIX])
   fi

   silo_LIBS="-L${silo_PREFIX}/lib ${silo_LIBS}"
fi

BTNG_AC_LOG(silo_INCLUDES is $silo_INCLUDES)
BTNG_AC_LOG(silo_LIBS is $silo_LIBS)

# END CASC_SUPPORT_SILO

])dnl End definition of CASC_SUPPORT_SILO
