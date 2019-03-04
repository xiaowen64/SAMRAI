#ifndef included_tbox_Raja_api
#define included_tbox_Raja_api

// #include "SAMRAI/pdat/ArrayData.h"

#include "RAJA/RAJA.hpp"

namespace SAMRAI {
namespace tbox {

namespace policy {
  struct sequential {};
  struct parallel {};
  struct host {};
}

namespace detail {

template <typename pol>
struct policy_traits{
  typedef void policy;
  typedef void raja_1d_policy;
  typedef void raja_2d_policy;
  typedef void raja_3d_policy;
};

template <>
struct policy_traits<policy::parallel> {
  typedef RAJA::cuda_exec_async<128> policy;

  using raja_1d_policy =
    RAJA::KernelPolicy<
          RAJA::statement::CudaKernelAsync<
            RAJA::statement::Tile<0, RAJA::statement::tile_fixed<128>, RAJA::cuda_block_x_loop,
              RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                RAJA::statement::Lambda<0>
              >
            >
          >
    >;

  using raja_2d_policy =
    RAJA::KernelPolicy<
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

  using raja_3d_policy =
    RAJA::KernelPolicy<
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

  typedef RAJA::cuda_reduce raja_reduction_policy;
};

template <>
struct policy_traits<policy::host> {
  typedef RAJA::seq_exec policy;

  typedef
    RAJA::KernelPolicy<
        RAJA::statement::For<0, RAJA::seq_exec, RAJA::statement::Lambda<0> >
    > raja_1d_policy;

  typedef
    RAJA::KernelPolicy<
        RAJA::statement::For<1, RAJA::seq_exec,
          RAJA::statement::For<0, RAJA::seq_exec, RAJA::statement::Lambda<0> > >
    > raja_2d_policy;
};

struct layout_traits {
   typedef RAJA::OffsetLayout<1, RAJA::Index_type> layout1d;
   typedef RAJA::OffsetLayout<2, RAJA::Index_type> layout2d;
   typedef RAJA::OffsetLayout<3, RAJA::Index_type> layout3d;
};

}

template <typename policy>
inline void
synchronize();

template<>
inline void
synchronize<tbox::policy::parallel>()
{
  RAJA::synchronize<RAJA::cuda_synchronize>();
}


/*
 * Simple forall
 */
template <typename policy, typename loop_body>
inline void for_all(int begin, int end, loop_body&& body) {
  RAJA::forall<typename detail::policy_traits<policy>::policy>(RAJA::RangeSegment(begin, end), body);
}



/*
 * for_all over hier::Box
 */
// template<typename policy, typename loop_body>
// inline void for_all(const hier::Box& box, loop_body body);

template<typename policy, typename loop_body>
inline void for_all1(const hier::Box& box, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

   RAJA::kernel<typename detail::policy_traits<policy>::raja_1d_policy >(
       RAJA::make_tuple(RAJA::RangeSegment(ifirst(0), ilast(0)+1)),
       body);
}

template<typename policy, typename loop_body>
inline void for_all1(const hier::Box& box, const int dim, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

   RAJA::kernel< typename detail::policy_traits<policy>::raja_1d_policy >(
       RAJA::make_tuple(RAJA::RangeSegment(ifirst(dim), ilast(dim)+1)),
       body);
}

template<typename policy, typename loop_body>
inline void for_all2(const hier::Box& box, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

   RAJA::kernel< typename detail::policy_traits<policy>::raja_2d_policy >(
       RAJA::make_tuple(RAJA::RangeSegment(ifirst(1), ilast(1)+1),
         RAJA::RangeSegment(ifirst(0), ilast(0)+1)),
       body);
}

template<typename policy, typename loop_body>
inline void for_all3(const hier::Box& box, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

  RAJA::kernel< typename detail::policy_traits<policy>::raja_3d_policy >(
      RAJA::make_tuple(RAJA::RangeSegment(ifirst(2), ilast(2)+1),
        RAJA::RangeSegment(ifirst(1), ilast(1)+1),
        RAJA::RangeSegment(ifirst(0), ilast(0)+1)),
      body);
}


template <size_t dim, typename T>
struct ArrayView {};

template <typename T>
struct ArrayView<1, T>  :
  public RAJA::View<T, detail::layout_traits::layout1d >
{
  using Layout = detail::layout_traits::layout1d;

  // T could be const?
  ArrayView<1, T>(T* data, const hier::Box& box, int depth = 0) :
     RAJA::View<T, Layout>(
         &data[depth * (box.size()-1)],
         RAJA::make_permuted_offset_layout(
           std::array<RAJA::Index_type, 1>{ {box.lower()[0]} }, std::array<RAJA::Index_type, 1>{ {box.upper()[0]} }, RAJA::as_array<RAJA::PERM_I>::get())){}

  ArrayView<1, T>(pdat::ArrayData<T>& data, int depth = 0) :
          ArrayView<1, T>(data.getPointer(depth), data.getBox()) {}
};

 template <typename T>
 struct ArrayView<2, T> : public RAJA::View<T, detail::layout_traits::layout2d >
 {
  using Layout = detail::layout_traits::layout2d;

   SAMRAI_INLINE ArrayView<2, T>(T* data, const hier::Box& box, int depth = 0) :
     RAJA::View<T, Layout>(
         &data[depth * (box.size())],
         RAJA::make_permuted_offset_layout(
           std::array<RAJA::Index_type, 2>{ {box.lower()[0], box.lower()[1]} },
           std::array<RAJA::Index_type, 2>{ {box.upper()[0], box.upper()[1]} },
           RAJA::as_array<RAJA::PERM_JI>::get())){}

   SAMRAI_INLINE ArrayView<2, T>(pdat::ArrayData<T>& data, int depth = 0) :
     ArrayView<2, T>(data.getPointer(depth), data.getBox()) {}
 };

 template <typename T>
 struct ArrayView<3, T> : public RAJA::View<T, detail::layout_traits::layout3d >
 {
  using Layout = detail::layout_traits::layout3d;

   SAMRAI_INLINE ArrayView<3, T>(T* data, const hier::Box& box, int depth = 0) :
     RAJA::View<T, Layout>(
         &data[depth * (box.size()-1)],
         RAJA::make_permuted_offset_layout(
           std::array<RAJA::Index_type, 3>{ {box.lower()[0], box.lower()[1], box.lower()[2]} },
           std::array<RAJA::Index_type, 3>{ {box.upper()[0], box.upper()[1], box.upper()[2]} },
           RAJA::as_array<RAJA::PERM_KJI>::get())){};

   SAMRAI_INLINE ArrayView<3, T>(pdat::ArrayData<T>& data, int depth = 0) :
     ArrayView<3, T>(data.getPointer(depth), data.getBox()) {}
 };

 template <typename T>
 struct ArrayView<1, const T>  :
   public RAJA::View<const T, detail::layout_traits::layout1d >
 {
   using Layout = detail::layout_traits::layout1d;

   // T could be const?
   ArrayView<1, const T>(const T* data, const hier::Box& box, int depth = 0) :
      RAJA::View<const T, Layout>(
          &data[depth * (box.size()-1)],
          RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 1>{ {box.lower()[0]} }, std::array<RAJA::Index_type, 1>{ {box.upper()[0]} }, RAJA::as_array<RAJA::PERM_I>::get())){}

   ArrayView<1, const T>(const pdat::ArrayData<T>& data, int depth = 0) :
     ArrayView<1, const T>(data.getPointer(depth), data.getBox()) {}
 };

template <typename T>
struct ArrayView<2, const T> :
  public RAJA::View<const T, detail::layout_traits::layout2d >
{
 using Layout = detail::layout_traits::layout2d;

  SAMRAI_INLINE ArrayView<2, const T>(const T* data, const hier::Box& box, int depth = 0) :
    RAJA::View<const T, Layout>(
        &data[depth * (box.size()-1)],
        RAJA::make_permuted_offset_layout(
          std::array<RAJA::Index_type, 2>{ {box.lower()[0], box.lower()[1]} },
          std::array<RAJA::Index_type, 2>{ {box.upper()[0], box.upper()[1]} },
          RAJA::as_array<RAJA::PERM_JI>::get())){}

  SAMRAI_INLINE ArrayView<2, const T>(const pdat::ArrayData<T>& data, int depth = 0) :
    ArrayView<2, const T>(data.getPointer(depth), data.getBox()) {}

};

template <typename T>
struct ArrayView<3, const T> :
  public RAJA::View<const T, detail::layout_traits::layout3d >
{
 using Layout = detail::layout_traits::layout3d;

  SAMRAI_INLINE ArrayView<3, const T>(const T* data, const hier::Box& box, int depth = 0) :
    RAJA::View<const T, Layout>(
        &data[depth * (box.size()-1)],
        RAJA::make_permuted_offset_layout(
          std::array<RAJA::Index_type, 3>{ {box.lower()[0], box.lower()[1], box.lower()[2]} },
          std::array<RAJA::Index_type, 3>{ {box.upper()[0], box.upper()[1], box.upper()[2]} },
          RAJA::as_array<RAJA::PERM_KJI>::get())){};

  SAMRAI_INLINE ArrayView<3, const T>(const pdat::ArrayData<T>& data, int depth = 0) :
    ArrayView<3, const T>(data.getPointer(depth), data.getBox()) {}
};

}
}

#endif
