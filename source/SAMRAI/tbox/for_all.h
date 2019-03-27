/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#ifndef included_tbox_for_all
#define included_tbox_for_all

#include "SAMRAI/SAMRAI_config.h"

#if defined(HAVE_RAJA)

#include "RAJA/RAJA.hpp"

#include <type_traits>
#include <tuple>
#include <cstdlib> // for std::size_t

namespace SAMRAI {
namespace tbox {

namespace policy {
struct sequential {};
struct parallel {};
}

namespace detail {

template <typename pol>
struct policy_traits {};

template <>
struct policy_traits<policy::sequential> {
   using Policy = RAJA::seq_exec;

   using Policy1d = RAJA::KernelPolicy<
      RAJA::statement::For<0, RAJA::seq_exec,
         RAJA::statement::Lambda<0>
      >
   >;

   using Policy2d = RAJA::KernelPolicy<
      RAJA::statement::For<1, RAJA::seq_exec,
         RAJA::statement::For<0, RAJA::seq_exec,
            RAJA::statement::Lambda<0>
         >
      >
   >;

   using Policy3d = RAJA::KernelPolicy<
      RAJA::statement::For<2, RAJA::seq_exec,
         RAJA::statement::For<1, RAJA::seq_exec,
            RAJA::statement::For<0, RAJA::seq_exec,
               RAJA::statement::Lambda<0>
            >
         >
      >
   >;

   using ReductionPolicy = RAJA::seq_reduce;
};

#if defined(HAVE_CUDA)

template <>
struct policy_traits<policy::parallel> {
   using Policy = RAJA::cuda_exec_async<128>;

   using Policy1d = RAJA::KernelPolicy<
      RAJA::statement::CudaKernelAsync<
         RAJA::statement::Tile<0, RAJA::statement::tile_fixed<128>, RAJA::cuda_block_x_loop,
            RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
               RAJA::statement::Lambda<0>
            >
         >
      >
   >;

   using Policy2d = RAJA::KernelPolicy<
      RAJA::statement::CudaKernelAsync<
         RAJA::statement::Tile<1, RAJA::statement::tile_fixed<32>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<0, RAJA::statement::tile_fixed<32>, RAJA::cuda_block_x_loop,
               RAJA::statement::For<1, RAJA::cuda_thread_y_loop,
                  RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                     RAJA::statement::Lambda<0>
                  >
               >
            >
         >
      >
   >;

   using Policy3d = RAJA::KernelPolicy<
      RAJA::statement::CudaKernelAsync<
         RAJA::statement::Tile<2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_z_loop,
            RAJA::statement::Tile<1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
               RAJA::statement::Tile<0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                  RAJA::statement::For<2, RAJA::cuda_thread_z_loop,
                     RAJA::statement::For<1, RAJA::cuda_thread_y_loop,
                        RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                           RAJA::statement::Lambda<0>
                        >
                     >
                  >
               >
            >
         >
      >
   >;

   using ReductionPolicy = RAJA::cuda_reduce;
};

#else

// TODO: Make this an OpenMP policy if that is defined
template <>
struct policy_traits<policy::parallel> {
   using Policy = RAJA::loop_exec;

   using Policy1d = RAJA::KernelPolicy<
      RAJA::statement::For<0, RAJA::loop_exec,
         RAJA::statement::Lambda<0>
      >
   >;

   using Policy2d = RAJA::KernelPolicy<
      RAJA::statement::For<1, RAJA::loop_exec,
         RAJA::statement::For<0, RAJA::loop_exec,
            RAJA::statement::Lambda<0>
         >
      >
   >;

   using Policy3d = RAJA::KernelPolicy<
      RAJA::statement::For<2, RAJA::loop_exec,
         RAJA::statement::For<1, RAJA::loop_exec,
            RAJA::statement::For<0, RAJA::loop_exec,
               RAJA::statement::Lambda<0>
            >
         >
      >
   >;

   using ReductionPolicy = RAJA::seq_reduce;
};

#endif // HAVE_CUDA

/*
tbox::parallel_for_all()  version that picks parallel policy (GPU if ENABLE_CUDA=ON)
tbox::for_all<pol>()      version that takes a custom RAJA policy (could be either host or device)
*/

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

template <>
struct for_all<1> {
   template <typename Policy, typename LoopBody>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename policy_traits<Policy>::Policy1d>(
         RAJA::make_tuple(make_range(ifirst, ilast, 0)),
         body);
   }
};

template <>
struct for_all<2> {
   template <typename Policy, typename LoopBody>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename policy_traits<Policy>::Policy2d>(
         RAJA::make_tuple(make_range(ifirst, ilast, 1),
                          make_range(ifirst, ilast, 0)),
         body);
   }
};

template <>
struct for_all<3> {
   template <typename Policy, typename LoopBody>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename policy_traits<Policy>::Policy3d>(
         RAJA::make_tuple(make_range(ifirst, ilast, 2),
                          make_range(ifirst, ilast, 1),
                          make_range(ifirst, ilast, 0)),
         body);
   }
};

} // namespace detail


template <typename Policy, typename LoopBody>
inline void for_all(int begin, int end, LoopBody body)
{
   RAJA::forall<typename detail::policy_traits<Policy>::Policy>(RAJA::RangeSegment(begin, end), body);
}

template <typename LoopBody>
inline void parallel_for_all(int begin, int end, LoopBody&& body)
{
   for_all<policy::parallel>(begin, end, body);
}


template<typename Policy, typename LoopBody>
inline void for_all(const hier::Box& box, const int dim, LoopBody body)
{
   RAJA::forall<typename detail::policy_traits<Policy>::Policy>(
      RAJA::RangeSegment(box.lower()(dim), box.upper()(dim)),
      body);
}

template<typename LoopBody>
inline void parallel_for_all(const hier::Box& box, const int dim, LoopBody body)
{
   for_all<policy::parallel>(box, dim, body);
}


template <typename Policy, typename LoopBody>
inline void for_all(const hier::Box& box, LoopBody body)
{
   constexpr int arg_count = detail::function_traits<LoopBody>::argument_count;
   detail::for_all<arg_count>::template eval<Policy>(box.lower(), box.upper(), body);
}

template <typename LoopBody>
inline void parallel_for_all(const hier::Box& box, LoopBody body)
{
   for_all<policy::parallel>(box, body);
}

// Reductions see https://raja.readthedocs.io/en/master/feature/reduction.html for these options
enum class Reduction {
   Sum,
   Min,
   Max,
   MinLoc,
   MaxLoc
};

// auto norm = tbox::parallel_reduction_variable<tbox::Reduction::Sum, double>::type
// or _t using trick
// auto norm = tbox::reduction_variable<Policy, tbox::Reduction::Sum, double>::type
// or _t using trick

template<typename Policy, Reduction R, typename TYPE = double>
struct reduction_variable;

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Sum, TYPE> {
   using type = RAJA::ReduceSum<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Min, TYPE> {
   using type = RAJA::ReduceMin<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Max, TYPE> {
   using type = RAJA::ReduceMax<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::MinLoc, TYPE> {
   using type = RAJA::ReduceMinLoc<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::MaxLoc, TYPE> {
   using type = RAJA::ReduceMaxLoc<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<Reduction R, typename TYPE = double>
using parallel_reduction_variable = reduction_variable<policy::parallel, R, TYPE>;

template<typename Policy, Reduction R, typename TYPE = double>
using reduction_variable_t = typename reduction_variable<Policy, R, TYPE>::type;

template<Reduction R, typename TYPE = double>
using parallel_reduction_variable_t = typename parallel_reduction_variable<R, TYPE>::type;

// Synchronization

template <typename policy>
inline void
synchronize() {}

template<>
inline void
synchronize<tbox::policy::parallel>()
{
   RAJA::synchronize<RAJA::cuda_synchronize>();
}

}
}

#endif

#endif // included_tbox_for_all
