#ifndef DEQN_AMR_RAJA_API_H
#define DEQN_AMR_RAJA_API_H

#include "SAMRAI/pdat/ArrayData.h"

#include "RAJA/RAJA.hxx"

namespace SAMRAI {
namespace tbox {

namespace policy {
  struct sequential {};
  struct parallel {};
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
  typedef RAJA::cuda_exec<128> policy;
  typedef RAJA::NestedPolicy< RAJA::ExecList< 
    RAJA::cuda_exec<128> > > raja_1d_policy;
  typedef RAJA::NestedPolicy< RAJA::ExecList< 
            RAJA::cuda_threadblock_y_exec<16>, 
            RAJA::cuda_threadblock_x_exec<16> > > raja_2d_policy;
  typedef RAJA::NestedPolicy< RAJA::ExecList<
            RAJA::cuda_threadblock_z_exec<8>, 
            RAJA::cuda_threadblock_y_exec<8>, 
            RAJA::cuda_threadblock_x_exec<8> > > raja_3d_policy;
};

struct layout_traits {
   typedef RAJA::OffsetLayout<1, int> layout1d;
   typedef RAJA::OffsetLayout<2, int> layout2d;
   typedef RAJA::OffsetLayout<3, int> layout3d;
};

}

/*
 * Simple forall
 */
template <typename policy, typename loop_body>
inline void for_all(int begin, int end, loop_body&& body) {
  RAJA::forall<typename detail::policy_traits<policy>::policy>(begin, end, body);
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

   RAJA::forallN< typename detail::policy_traits<policy>::raja_1d_policy >(
       RAJA::RangeSegment(ifirst(0), ilast(0)+1),
       body);
}

template<typename policy, typename loop_body>
inline void for_all2(const hier::Box& box, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

   RAJA::forallN< typename detail::policy_traits<policy>::raja_2d_policy >(
       RAJA::RangeSegment(ifirst(1), ilast(1)+1),
       RAJA::RangeSegment(ifirst(0), ilast(0)+1), 
       body);
}

template<typename policy, typename loop_body>
inline void for_all3(const hier::Box& box, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

  RAJA::forallN< typename detail::policy_traits<policy>::raja_3d_policy >(
      RAJA::RangeSegment(ifirst(2), ilast(2)+1),
      RAJA::RangeSegment(ifirst(1), ilast(1)+1),
      RAJA::RangeSegment(ifirst(0), ilast(0)+1), 
      body);
}


template <size_t dim, typename T>
struct ArrayView {};

template <typename T>
struct ArrayView<1, T>  : 
  public RAJA::View<T, detail::layout_traits::layout1d >
{
  using Layout = detail::layout_traits::layout1d;

  ArrayView<1, T>(pdat::ArrayData<T>& data, int depth = 0) :
    RAJA::View<T, Layout >(
        data.getPointer(depth),
        RAJA::make_permuted_offset_layout(
          {data.getBox().lower()[0]}, {data.getBox().upper()[0]}, RAJA::PERM_I::value)){}

  // T could be const?
  ArrayView<1, T>(T* data, const hier::Box& box, int depth = 0) :
     RAJA::View<T, Layout>(
         &data[depth * (box.size()-1)],
         RAJA::make_permuted_offset_layout(
           {box.lower()[0]}, {box.upper()[0]}, RAJA::PERM_I::value)){}
};

 template <typename T>
 struct ArrayView<2, T> : public RAJA::View<T, detail::layout_traits::layout2d >
 {
  using Layout = detail::layout_traits::layout2d;

   SAMRAI_INLINE ArrayView<2, T>(pdat::ArrayData<T>& data, int depth = 0) :
     RAJA::View<T, Layout>(
         data.getPointer(depth),
         RAJA::make_permuted_offset_layout(
           {data.getBox().lower()[0], data.getBox().lower()[1]}, 
           {data.getBox().upper()[0], data.getBox().upper()[1]},
           RAJA::PERM_JI::value)){}

   SAMRAI_INLINE ArrayView<2, T>(T* data, const hier::Box& box, int depth = 0) :
     RAJA::View<T, Layout>(
         &data[depth * (box.size()-1)],
         RAJA::make_permuted_offset_layout(
           {box.lower()[0], box.lower()[1]}, 
           {box.upper()[0], box.upper()[1]},
           RAJA::PERM_JI::value)){}
 };
 
 template <typename T>
 struct ArrayView<3, T> : public RAJA::View<T, detail::layout_traits::layout3d >
 {
  using Layout = detail::layout_traits::layout3d;

   SAMRAI_INLINE ArrayView<3, T>(pdat::ArrayData<T>& data, int depth = 0) :
     RAJA::View<T, Layout>(
         data.getPointer(depth),
         RAJA::make_permuted_offset_layout(
           {data.getBox().lower()[0], data.getBox().lower()[1], data.getBox().lower()[2]}, 
           {data.getBox().upper()[0], data.getBox().upper()[1], data.getBox().upper()[2]},
           RAJA::PERM_KJI::value)) {}

   SAMRAI_INLINE ArrayView<3, T>(T* data, const hier::Box& box, int depth = 0) :
     RAJA::View<T, Layout>(
         &data[depth * (box.size()-1)],
         RAJA::make_permuted_offset_layout(
           {box.lower()[0], box.lower()[1], box.lower()[2]}, 
           {box.upper()[0], box.upper()[1], box.upper()[2]},
           RAJA::PERM_KJI::value)){};
 };

 template <typename T>
 struct ArrayView<1, const T>  : 
   public RAJA::View<const T, detail::layout_traits::layout1d >
 {
   using Layout = detail::layout_traits::layout1d;
 
   ArrayView<1, const T>(const pdat::ArrayData<T>& data, int depth = 0) :
     RAJA::View<const T, Layout >(
         data.getPointer(depth),
         RAJA::make_permuted_offset_layout(
           {data.getBox().lower()[0]}, {data.getBox().upper()[0]}, RAJA::PERM_I::value)){}
 
   // T could be const?
   ArrayView<1, const T>(const T* data, const hier::Box& box, int depth = 0) :
      RAJA::View<const T, Layout>(
          &data[depth * (box.size()-1)],
          RAJA::make_permuted_offset_layout(
            {box.lower()[0]}, {box.upper()[0]}, RAJA::PERM_I::value)){}
 };

template <typename T>
struct ArrayView<2, const T> : 
  public RAJA::View<const T, detail::layout_traits::layout2d >
{
 using Layout = detail::layout_traits::layout2d;

  SAMRAI_INLINE ArrayView<2, const T>(const pdat::ArrayData<T>& data, int depth = 0) :
    RAJA::View<const T, Layout>(
        data.getPointer(depth),
        RAJA::make_permuted_offset_layout(
          {data.getBox().lower()[0], data.getBox().lower()[1]}, 
          {data.getBox().upper()[0], data.getBox().upper()[1]},
          RAJA::PERM_JI::value)){}

  SAMRAI_INLINE ArrayView<2, const T>(const T* data, const hier::Box& box, int depth = 0) :
    RAJA::View<const T, Layout>(
        &data[depth * (box.size()-1)],
        RAJA::make_permuted_offset_layout(
          {box.lower()[0], box.lower()[1]}, 
          {box.upper()[0], box.upper()[1]},
          RAJA::PERM_JI::value)){}
};

template <typename T>
struct ArrayView<3, const T> : 
  public RAJA::View<const T, detail::layout_traits::layout3d >
{
 using Layout = detail::layout_traits::layout3d;

  SAMRAI_INLINE ArrayView<3, const T>(const pdat::ArrayData<T>& data, int depth = 0) :
    RAJA::View<const T, Layout>(
        data.getPointer(depth),
        RAJA::make_permuted_offset_layout(
          {data.getBox().lower()[0], data.getBox().lower()[1], data.getBox().lower()[2]}, 
          {data.getBox().upper()[0], data.getBox().upper()[1], data.getBox().upper()[2]},
          RAJA::PERM_KJI::value)) {}

  SAMRAI_INLINE ArrayView<3, const T>(const T* data, const hier::Box& box, int depth = 0) :
    RAJA::View<const T, Layout>(
        &data[depth * (box.size()-1)],
        RAJA::make_permuted_offset_layout(
          {box.lower()[0], box.lower()[1], box.lower()[2]}, 
          {box.upper()[0], box.upper()[1], box.upper()[2]},
          RAJA::PERM_KJI::value)){};
};
 
 // template <typename T, int DIM>
 // struct CellView : public ArrayView<T, DIM>
 // {
 //   SAMRAI_INLINE CellView(const boost::shared_ptr<pdat::CellData<T> >& data)
 //     : ArrayView<
 //   {
 //   }
 // };

}
}

#endif //DEQN_AMR_RAJA_API_H
