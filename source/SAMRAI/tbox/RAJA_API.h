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

template <typename policy>
struct policy_traits{
  typedef void policy;
  typedef void raja_1d_policy;
  typedef void raja_2d_policy;
  typedef void raja_3d_policy;
};

struct policy_traits<policy::parallel> {
  typedef RAJA::cuda_exec<128> policy;
  typedef RAJA::NestedPolicy< RAJA::ExecList< 
    RAJA::cuda_exec<256> > > raja_1d_policy;
  typedef RAJA::NestedPolicy< RAJA::ExecList< 
            RAJA::cuda_exec<256>, 
            RAJA::cuda_exec<256> > > raja_2d_policy
  typedef RAJA::NestedPolicy< RAJA::ExecList<
            RAJA::cuda_exec<256>, 
            RAJA::cuda_exec<256>, 
            RAJA::cuda_exec<256> > > raja_3d_policy;
};

/*
 * Simple forall
 */
template <typename policy, typename loop_body>
inline void for_all(int begin, int end, loop_body&& body) {
  RAJA::forall<detail::policy_traits<policy>::policy>(begin, end, body);
}


/*
 * for_all over hier::Box
 */
template<typename policy, int dimension, typename loop_body>
void for_all(const hier::Box& box, loop_body body);

template<typename policy, typename loop_body>
inline void for_all<1>(const hier::Box& box, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

   RAJA::forallN< detail::policy_traits<policy>::raja_1d_policy >(
       RAJA::RangeSegment(ifirst(1), ilast(1)+1), 
       RAJA::RangeSegment(ifirst(0), ilast(0)+1),
       body);
}

template<typename policy, typename loop_body>
inline void for_all<2>(const hier::Box& box, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

   RAJA::forallN< detail::policy_traits<policy>::raja_2d_policy >(
       RAJA::RangeSegment(ifirst(1), ilast(1)+1),
       RAJA::RangeSegment(ifirst(0), ilast(0)+1), 
       body);
}

template<typename policy, typename loop_body>
inline void for_all<3>(const hier::Box& box, loop_body body)
{
  const hier::Index ifirst = box.lower();
  const hier::Index ilast = box.upper();

  RAJA::forallN< detail::policy_traits<policy>::raja_3d_policy >(
      RAJA::RangeSegment(ifirst(2), ilast(2)+1),
      RAJA::RangeSegment(ifirst(1), ilast(1)+1),
      RAJA::RangeSegment(ifirst(0), ilast(0)+1), 
      body);
}

template <int DIM>
struct Layout {};

struct Layout<1> : RAJA::OffsetLayout<int, RAJA::PERM_I, int> {};

struct Layout<2> : RAJA::OffsetLayout<int, RAJA::PERM_JI, int, int> {};

struct Layout<3> : RAJA::OffsetLayout<int, RAJA::PERM_KJI, int, int, int> {};

template <typename T, int DIM>
struct ArrayView : public RAJA::View<T, Layout<DIM> >

template <typename T>
struct ArrayView<1> : public RAJA::View<T, Layout<1> >
{
  SAMRAI_INLINE ArrayView(const ArrayData<T>& data, int depth = 0) :
    RAJA::View<T, Layout>(
        data.getPointer(depth),
        Layout({data.getBox().lower()[0]}, {data.getBox().upper()[0]})
});

template <typename T>
struct ArrayView<2> : public RAJA::View<T, Layout<2> >
{
  SAMRAI_INLINE ArrayView(const ArrayData& data, int depth = 0) :
    RAJA::View<T, Layout>(
        data.getPointer(depth),
        Layout({data.getBox().lower()[0], data.getBox().lower()[1]}, 
          {data.getBox().upper()[0], data.getBox().upper()[1]}) {}
});

template <typename T>
struct ArrayView<3> : public RAJA::View<T, Layout<3> >
{
  SAMRAI_INLINE ArrayView(const ArrayData& data, int depth = 0) :
    RAJA::View<T, Layout>(
        data.getPointer(depth),
        Layout({data.getBox().lower()[0], data.getBox().lower()[1], data.getBox().lower()[2]}, 
          {data.getBox().upper()[0], data.getBox().upper()[1], data.getBox().upper()[2]}})
        {
        }
});

/* 
 * ArrayView constructed with a buffer and an exsting Box, used when packing
 * data into a buffer.
 */
template <typename T>
struct ArrayView<1> : public RAJA::View<T, Layout<1> >
{
  SAMRAI_INLINE ArrayView(T* data, const hier::Box& box, int depth = 0) :
    RAJA::View<T, Layout>(
        &data[depth * (box.size()-1)],
        Layout({box.lower()[0]}, {box.upper()[0]})
});

template <typename T>
struct ArrayView<2> : public RAJA::View<T, Layout<2> >
{
  SAMRAI_INLINE ArrayView(T* data, const hier::Box& box, int depth = 0) :
    RAJA::View<T, Layout>(
        box[depth * (box.size()-1)],
        Layout({box.lower()[0], box.lower()[1]}, 
          {box.upper()[0], box.upper()[1]}) {}
});

template <typename T>
struct ArrayView<3> : public RAJA::View<T, Layout<3> >
{
  SAMRAI_INLINE ArrayView(T* data, const hier::Box& box, int depth = 0) :
    RAJA::View<T, Layout>(
        &data[depth * (box.size()-1)],
        Layout({box.lower()[0], box.lower()[1], box.lower()[2]}, 
          {box.upper()[0], box.upper()[1], box.upper()[2]}})
        {
        }
});

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
