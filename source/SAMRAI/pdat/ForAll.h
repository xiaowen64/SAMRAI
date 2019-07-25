/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#ifndef included_pdat_ForAll
#define included_pdat_ForAll

#include "SAMRAI/SAMRAI_config.h"

#if defined(HAVE_RAJA)
#include "RAJA/RAJA.hpp"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/ExecutionPolicy.h"

#include <type_traits>
#include <tuple>
#include <cstdlib> // for std::size_t

namespace SAMRAI {
namespace pdat {

/*
parallel_for_all()  version that picks parallel policy (GPU if ENABLE_CUDA=ON)
for_all<policy>()      version that takes a custom RAJA policy (could be either host or device)
*/

namespace detail {

template<typename T>
struct function_traits : function_traits<decltype(&T::operator())> {};

// free function
template<typename R, typename... Args>
struct function_traits<R(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// pointer to function
template<typename R, typename... Args>
struct function_traits<R (*)(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// member function
template<typename T, typename R, typename... Args>
struct function_traits<R (T::*)(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// const member function
template<typename T, typename R, typename... Args>
struct function_traits<R (T::*)(Args...) const> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

inline RAJA::RangeSegment make_range(const hier::Index& ifirst, const hier::Index& ilast, std::size_t index)
{ return RAJA::RangeSegment(ifirst(index), ilast(index)+1); }

template <int ArgumentCount>
struct for_all {};

template <int ArgumentCount>
struct for_all_x {};

template <>
struct for_all<1> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy1d>(
         RAJA::make_tuple(make_range(ifirst, ilast, 0)),
         body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
         RAJA::make_tuple(make_range(ifirst, ilast, 0)),
         body);
   }
};

template <>
struct for_all<2> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy2d>(
         RAJA::make_tuple(make_range(ifirst, ilast, 1),
                          make_range(ifirst, ilast, 0)),
         body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
         RAJA::make_tuple(make_range(ifirst, ilast, 1),
                          make_range(ifirst, ilast, 0)),
         body);
   }
};

template <>
struct for_all_x<2> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy2d>(
         RAJA::make_tuple(make_range(ifirst, ilast, 0),
                          make_range(ifirst, ilast, 1)),
         body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
         RAJA::make_tuple(make_range(ifirst, ilast, 0),
                          make_range(ifirst, ilast, 1)),
         body);
   }
};

template <>
struct for_all<3> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy3d>(
         RAJA::make_tuple(make_range(ifirst, ilast, 2),
                          make_range(ifirst, ilast, 1),
                          make_range(ifirst, ilast, 0)),
         body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
         RAJA::make_tuple(make_range(ifirst, ilast, 2),
                          make_range(ifirst, ilast, 1),
                          make_range(ifirst, ilast, 0)),
         body);
   }
};

template <>
struct for_all_x<3> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy3d>(
         RAJA::make_tuple(make_range(ifirst, ilast, 0),
                          make_range(ifirst, ilast, 1),
                          make_range(ifirst, ilast, 2)),
         body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
         RAJA::make_tuple(make_range(ifirst, ilast, 2),
                          make_range(ifirst, ilast, 1),
                          make_range(ifirst, ilast, 0)),
         body);
   }
};

} // namespace detail

// does NOT include end
template <typename Policy, typename LoopBody,
          typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
inline void for_all(int begin, int end, LoopBody body)
{
   RAJA::forall<typename tbox::detail::policy_traits<Policy>::Policy>(RAJA::RangeSegment(begin, end), body);
}

template <typename Policy, typename LoopBody,
          typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
inline void for_all(int begin, int end, LoopBody body)
{
   RAJA::forall<Policy>(RAJA::RangeSegment(begin, end), body);
}

// does NOT include end
template <typename LoopBody>
inline void parallel_for_all(int begin, int end, LoopBody body)
{
   for_all<tbox::policy::parallel>(begin, end, body);
}


template<typename Policy, typename LoopBody>
inline void for_all(const hier::Box& box, const int dim, LoopBody body)
{
   for_all<Policy>(box.lower()(dim), box.upper()(dim)+1, body);
}

template<typename LoopBody>
inline void parallel_for_all(const hier::Box& box, const int dim, LoopBody body)
{
   for_all<tbox::policy::parallel>(box.lower()(dim), box.upper()(dim)+1, body);
}


template <typename Policy, typename LoopBody>
inline void for_all(const hier::Box& box, LoopBody body)
{
   constexpr int arg_count = detail::function_traits<LoopBody>::argument_count;
   detail::for_all<arg_count>::template eval<Policy>(box.lower(), box.upper(), body);
}

template <typename Policy, typename LoopBody>
inline void for_all_x(const hier::Box& box, LoopBody body)
{
   constexpr int arg_count = detail::function_traits<LoopBody>::argument_count;
   detail::for_all_x<arg_count>::template eval<Policy>(box.lower(), box.upper(), body);
}

template <typename LoopBody>
inline void parallel_for_all(const hier::Box& box, LoopBody body)
{
   for_all<tbox::policy::parallel>(box, body);
}

template <typename LoopBody>
inline void parallel_for_all_x(const hier::Box& box, LoopBody body)
{
   for_all_x<tbox::policy::parallel>(box, body);
}

}
}

#endif

#endif // included_pdat_ForAll
