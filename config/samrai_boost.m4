dnl Define a macro for supporting BOOST

AC_DEFUN([CASC_SUPPORT_BOOST],[

# Begin CASC_SUPPORT_BOOST
# Defines boost_PREFIX boost_INCLUDES and boost_LIBS.
AC_ARG_WITH(boost,
[ --with-boost[=PATH]  Use BOOST and specify where BOOST is installed.],
, with_boost=no)

case "$with_boost" in
  no)
    AC_MSG_NOTICE([configuring without BOOST support])
    : Do nothing
  ;;
  yes)
    # BOOST install path was not specified.
    # Look in a couple of standard locations to probe if 
    # BOOST header files are there.
    AC_MSG_CHECKING([for BOOST installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/boost/tr1/unordered_map.hpp; then
        boost_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$boost_PREFIX])
  ;;
  *)
    # BOOST install path was specified.
    AC_MSG_CHECKING([for BOOST installation])

    if test -f ${with_boost}/include/boost/tr1/unordered_map.hpp; then
        boost_PREFIX=$with_boost
        boost_INCLUDES="-I${boost_PREFIX}/include"
        boost_LIBS="-L${boost_PREFIX}/lib -lboost_regex"
        AC_MSG_RESULT([$boost_PREFIX])
    else
        AC_MSG_RESULT([$boost_PREFIX])
        AC_MSG_ERROR([BOOST not found in $with_boost])
    fi
  ;;
esac



# Test compiling an BOOST application

# NOTE that AC_SEARCH_LIBS didn't work completely so
# use a more complicated example program to see
# if that will catch when BOOST is not working.
if test "${boost_PREFIX+set}" = set; then

   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether BOOST link works)
   AC_LANG_PUSH(C++)
   CASC_PUSH_COMPILER_STATE
   # NOTE lib z and m were from BTNG macro.
   LIBS="${LIBS} ${boost_LIBS}"
   CXXFLAGS="${CXXFLAGS} ${boost_INCLUDES}"
   AC_LINK_IFELSE([
      #include <boost/tr1/unordered_map.hpp>
      #include <boost/tr1/regex.hpp>
      #include <string>

      int main() {
         boost::unordered_map<std::string, int> my_map;

         my_map[["one"]] = 1;
         my_map[["two"]] = 2;
         std::string expr1 = ["samrai"];
         const boost::regex e(expr1);
         return 0;
      }
      ], 
      casc_boost_compile=yes,
      casc_boost_compile=no)
   CASC_POP_COMPILER_STATE
   AC_LANG_POP
   AC_MSG_RESULT($casc_boost_compile)

   if test "$casc_boost_compile" = no; then
      AC_MSG_ERROR([BOOST compile/link test failed])
   fi
fi

# END CASC_SUPPORT_BOOST

])dnl End definition of CASC_SUPPORT_BOOST

