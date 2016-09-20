include(FortranCInterface)
FortranCInterface_HEADER(
  ${CMAKE_BINARY_DIR}/include/SAMRAI/FC.h
  MACRO_NAMESPACE "CMAKE_FORTRAN_")

include(CheckIncludeFiles)

# /* Define if bool type is not properly supported */
# #undef BOOL_IS_BROKEN
# 
if (${ENABLE_BOX_COUNTING})
  set (BOX_TELEMETRY On)
endif()
# 
# /* Define if std::fill_n returns void */
# #undef CASC_STD_FILL_N_RETURNS_VOID
# 
# /* Define if DBL_MAX is not in float.h */
# #undef DBL_MAX_IS_BROKEN
# 
# /* Define if DBL_SNAN is not in float.h */
# #undef DBL_SNAN_IS_BROKEN
# 
# /* Enable assertion checking */
# #undef DEBUG_CHECK_ASSERTIONS
# 
# /* Enable SAMRAI developer assertion checking */
# #undef DEBUG_CHECK_DEV_ASSERTIONS
# 
# /* Enable assertion checking for dimensions */
# #undef DEBUG_CHECK_DIM_ASSERTIONS
# 
# /* Initialize new memory to undefined values in debug mode */
# #undef DEBUG_INITIALIZE_UNDEFINED
# 
if (${ENABLE_TIMERS})
  set(ENABLE_SAMRAI_TIMERS On)
endif ()

# 
# /* Define to dummy `main' function (if any) required to link to the Fortran
#    libraries. */
# #undef F77_DUMMY_MAIN
# 
# /* Define to a macro mangling the given C identifier (in lower and upper
#    case), which must not contain underscores, for linking with Fortran. */
# #define SAMRAI_F77_FUNC CMAKE_FORTRAN_GLOBAL
# 
# /* As SAMRAI_F77_FUNC, but for C identifiers containing underscores. */
# #define SAMRAI_F77_FUNC_ CMAKE_FORTRAN_GLOBAL_
# 
# /* Define if F77 and FC dummy `main' functions are identical. */
# #undef FC_DUMMY_MAIN_EQ_F77
# 
# /* Define if FLT_MAX is not in float.h */
# #undef FLT_MAX_IS_BROKEN
# 
# /* Define if FLT_SNAN is not in float.h */
# #undef FLT_SNAN_IS_BROKEN
# 
# /* BLAS library is available so use it */
# #undef HAVE_BLAS
# 
# /* BOOST headers are available to use */
# #undef HAVE_BOOST_HEADERS
# 
# /* HAVE_CMATH */
# #undef HAVE_CMATH
# 
# /* HAVE_CMATH_ISNAN */
# #undef HAVE_CMATH_ISNAN
# 
# /* HAVE_CTIME */
# #undef HAVE_CTIME
# 
# /* HAVE_EXCEPTION_HANDLING */
# #undef HAVE_EXCEPTION_HANDLING
# 
# /* HDF5 library is available so use it */
# #undef HAVE_HDF5
# 
# /* HYPRE library is available so use it */
# #undef HAVE_HYPRE
# 
# /* HAVE_INLINE_ISNAND */
# #undef HAVE_INLINE_ISNAND
# 
# /* Define to 1 if you have the <inttypes.h> header file. */
# #undef HAVE_INTTYPES_H
# 
# /* HAVE_IOMANIP_LEFT */
# #undef HAVE_IOMANIP_LEFT
# 
# /* HAVE_ISNAN */
# #undef HAVE_ISNAN
# 
# /* HAVE_ISNAND */
# #undef HAVE_ISNAND
# 
# /* HAVE_ISNAN_TEMPLATE */
# #undef HAVE_ISNAN_TEMPLATE
# 
# /* HAVE_ISO_SSTREAM */
# #undef HAVE_ISO_SSTREAM
# 
# /* LAPACK library is available so use it */
# #undef HAVE_LAPACK
# 
# /* Define to 1 if you have the `z' library (-lz). */
# #undef HAVE_LIBZ
# 
# /* Define if you have the 'mallinfo' function. */
# #undef HAVE_MALLINFO
# 
# check_include_files(malloc.h HAVE_MALLOC_H)
# 
# /* HAVE_MEMBER_FUNCTION_SPECIALIZATION */
# #undef HAVE_MEMBER_FUNCTION_SPECIALIZATION

check_include_files(memory.h HAVE_MEMORY_H)

# /* MPI library is present */
# #undef HAVE_MPI
# 
# /* HAVE_NAMESPACE */
# #undef HAVE_NAMESPACE
# 
# /* HAVE_NEW_PLACEMENT_OPERATOR */
# #undef HAVE_NEW_PLACEMENT_OPERATOR
# 
# /* OPENMP is available */
# #undef HAVE_OPENMP
# 
# /* PETSC library is available so use it */
# #undef HAVE_PETSC
# 
# /* HAVE_PRAGMA_STATIC_DATA_SPECIALIZATION */
# #undef HAVE_PRAGMA_STATIC_DATA_SPECIALIZATION
# 
# /* PTSCOTCH headers are available to use */
# #undef HAVE_PTSCOTCH
# 
# /* SILO library is available so use it */
# #undef HAVE_SILO
# 
# /* HAVE_SSTREAM */
# #undef HAVE_SSTREAM
# 
# /* HAVE_STANDARD_STATIC_DATA_SPECIALIZATION */
# #undef HAVE_STANDARD_STATIC_DATA_SPECIALIZATION
# 
# /* HAVE_STATIC_DATA_INSTANTIATION */
# #undef HAVE_STATIC_DATA_INSTANTIATION

check_include_files(stdint.h HAVE_STDINT_H)

check_include_files(stdlib.h HAVE_STDLIB_H)

# /* Define to 1 if cpp supports the ANSI # stringizing operator. */
# #undef HAVE_STRINGIZE

check_include_files(strings.h HAVE_STRINGS_H)

check_include_files(string.h HAVE_STRING_H)

# /* HAVE_SUNDIALS */
# #undef HAVE_SUNDIALS

check_include_files(sys/stat.h HAVE_SYS_STAT_H)

check_include_files(sys/times.h HAVE_SYS_TIMES_H)

check_include_files(sys/types.h HAVE_SYS_TYPES_H)

# /* HAVE_TAU */
# #undef HAVE_TAU
# 
# /* Thread Building Blocks are available to use */
# #undef HAVE_TBB
# 
# /* HAVE_TEMPLATE_COMPLEX */
# #undef HAVE_TEMPLATE_COMPLEX

check_include_files(unistd.h HAVE_UNISTD_H)

# /* HAVE_VAMPIR */
# #undef HAVE_VAMPIR
# 
# /* X11 library is present */
# #undef HAVE_X11
# 
# /* "Compiling with XDR support" */
# #undef HAVE_XDR
# 
# /* Define if the host system is Solaris */
# #undef HOST_OS_IS_SOLARIS
# 
# /* Hypre library is configured for sequential mode */
# #undef HYPRE_SEQUENTIAL
# 
if (${ENABLE_DEPRECATED})
  set (INCLUDE_DEPRECATED On)
endif ()
# 
# /* Header file for iomanip */
# #define IOMANIP_HEADER_FILE <iomanip>
# 
# /* The iomanip header file is broken */
# #undef IOMANIP_IS_BROKEN
# 
# /* Header file for iostream */
# #undef IOSTREAM_HEADER_FILE
# 
# /* The iostream header file is broken */
# #undef IOSTREAM_IS_BROKEN
# 
# /* LACKS_CMATH */
# #undef LACKS_CMATH
# 
# /* LACKS_CMATH_ISNAN */
# #undef LACKS_CMATH_ISNAN
# 
# /* LACKS_CTIME */
# #undef LACKS_CTIME
# 
# /* LACKS_EXCEPTION_HANDLING */
# #undef LACKS_EXCEPTION_HANDLING
# 
# /* Hypre library is missing */
# #undef LACKS_HYPRE
# 
# /* LACKS_INLINE_ISNAND */
# #undef LACKS_INLINE_ISNAND
# 
# /* LACKS_IOMANIP_LEFT */
# #undef LACKS_IOMANIP_LEFT
# 
# /* LACKS_ISNAN */
# #undef LACKS_ISNAN
# 
# /* LACKS_ISNAND */
# #undef LACKS_ISNAND
# 
# /* LACKS_ISNAN_TEMPLATE */
# #undef LACKS_ISNAN_TEMPLATE
# 
# /* LACKS_MEMBER_FUNCTION_SPECIALIZATION */
# #undef LACKS_MEMBER_FUNCTION_SPECIALIZATION
# 
# /* MPI library is missing */
# #define LACKS_MPI
# 
# /* LACKS_NAMESPACE */
# #undef LACKS_NAMESPACE
# 
# /* LACKS_NEW_PLACEMENT_OPERATOR */
# #undef LACKS_NEW_PLACEMENT_OPERATOR
# 
# /* LACKS_PRAGMA_STATIC_DATA_SPECIALIZATION */
# #undef LACKS_PRAGMA_STATIC_DATA_SPECIALIZATION
# 
# /* LACKS_PROPER_XDR_HEADER */
# #undef LACKS_PROPER_XDR_HEADER
# 
# /* LACKS_SSTREAM */
# #undef LACKS_SSTREAM
# 
# /* LACKS_STANDARD_STATIC_DATA_SPECIALIZATION */
# #undef LACKS_STANDARD_STATIC_DATA_SPECIALIZATION
# 
# /* LACKS_STATIC_DATA_INSTANTIATION */
# #undef LACKS_STATIC_DATA_INSTANTIATION
# 
# /* LACKS_SUNDIALS */
# #undef LACKS_SUNDIALS
# 
# /* LACKS_TAU */
# #undef LACKS_TAU
# 
# /* LACKS_TEMPLATE_COMPLEX */
# #undef LACKS_TEMPLATE_COMPLEX
# 
# /* LACKS_VAMPIR */
# #undef LACKS_VAMPIR
# 
# /* X11 library is missing */
# #undef LACKS_X11
# 
# /* LACK_ISO_SSTREAM */
# #undef LACK_ISO_SSTREAM
# 
# /* Define if NAN is not in float.h */
# #undef NAN_IS_BROKEN
# 
# /* Optimized build */
# #undef OPT_BUILD
# 
# /* The type ostringstream is broken */
# #undef OSTRINGSTREAM_TYPE_IS_BROKEN
# 
# /* The type ostrstream is broken */
# #define OSTRSTREAM_TYPE_IS_BROKEN
# 
# /* Define if restrict is not properly supported */
# #undef RESTRICT_IS_BROKEN

set(SAMRAI_MAXIMUM_DIMENSION ${MAXDIM})

# 
# /* Define to 1 if you have the ANSI C header files. */
# #undef STDC_HEADERS
# 
# /* Header file for stl-sstream */
# #undef STL_SSTREAM_HEADER_FILE
# 
# /* The stl-sstream header file is broken */
# #undef STL_SSTREAM_IS_BROKEN
# 
# /* Define to 1 if the X Window System is missing or not being used. */
# #undef X_DISPLAY_MISSING
# 
# /* Kludgey thing inserted by configure.in */
# #undef _POWER
# 
# /* Configure for compiling on BGL family of machines */
# #undef __BGL_FAMILY__
# 
# /*
#  * Prevent inclusion of mpi C++ bindings in mpi.h includes.
#  * This is done in here rather than SAMRAI_MPI.h since other
#  * files include MPI.h, such as PETSc and hypre.
#  */
# #ifndef MPI_NO_CPPBIND
# #define MPI_NO_CPPBIND
# #endif
# 
# #ifndef MPICH_SKIP_MPICXX
# #define MPICH_SKIP_MPICXX
# #endif
# 
# #ifndef OMPI_SKIP_MPICXX
# #define OMPI_SKIP_MPICXX
# #endif

configure_file(${PROJECT_SOURCE_DIR}/config/SAMRAI_config.h.cmake.in ${CMAKE_BINARY_DIR}/include/SAMRAI/SAMRAI_config.h)
