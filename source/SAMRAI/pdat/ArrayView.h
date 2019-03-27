/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Templated array data structure supporting patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayView
#define included_pdat_ArrayView

#include "SAMRAI/SAMRAI_config.h"

#if defined(HAVE_RAJA)

#include "RAJA/RAJA.hpp"

namespace SAMRAI {
namespace pdat {

namespace detail {

struct layout_traits {
   using Layout1d = RAJA::OffsetLayout<1, RAJA::Index_type>;
   using Layout2d = RAJA::OffsetLayout<2, RAJA::Index_type>;
   using Layout3d = RAJA::OffsetLayout<3, RAJA::Index_type>;
};

} // namespace detail

template<int DIM, class TYPE>
struct ArrayView {};

template<class TYPE>
struct ArrayView<1, TYPE> : public RAJA::View<TYPE, detail::layout_traits::Layout1d>
{
   using Layout = detail::layout_traits::Layout1d;

   ArrayView<1, TYPE>(TYPE* data, const hier::Box& box, int depth = 0) :
      RAJA::View<TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 1>{ {box.lower()[0]} },
            std::array<RAJA::Index_type, 1>{ {box.upper()[0]} },
            RAJA::as_array<RAJA::PERM_I>::get())){}
};

template<class TYPE>
struct ArrayView<2, TYPE> : public RAJA::View<TYPE, detail::layout_traits::Layout2d>
{
   using Layout = detail::layout_traits::Layout2d;

   SAMRAI_INLINE ArrayView<2, TYPE>(TYPE* data, const hier::Box& box, int depth = 0) :
      RAJA::View<TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 2>{ {box.lower()[0], box.lower()[1]} },
            std::array<RAJA::Index_type, 2>{ {box.upper()[0], box.upper()[1]} },
            RAJA::as_array<RAJA::PERM_JI>::get())){}
};

template<class TYPE>
struct ArrayView<3, TYPE> : public RAJA::View<TYPE, detail::layout_traits::Layout3d>
{
   using Layout = detail::layout_traits::Layout3d;

   SAMRAI_INLINE ArrayView<3, TYPE>(TYPE* data, const hier::Box& box) :
      RAJA::View<TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 3>{ {box.lower()[0], box.lower()[1], box.lower()[2]} },
            std::array<RAJA::Index_type, 3>{ {box.upper()[0], box.upper()[1], box.upper()[2]} },
            RAJA::as_array<RAJA::PERM_KJI>::get())){};
};

template<class TYPE>
struct ArrayView<1, const TYPE> : public RAJA::View<const TYPE, detail::layout_traits::Layout1d>
{
   using Layout = detail::layout_traits::Layout1d;

   ArrayView<1, const TYPE>(const TYPE* data, const hier::Box& box) :
      RAJA::View<const TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 1>{ {box.lower()[0]} },
            std::array<RAJA::Index_type, 1>{ {box.upper()[0]} },
            RAJA::as_array<RAJA::PERM_I>::get())){}
};

template<class TYPE>
struct ArrayView<2, const TYPE> : public RAJA::View<const TYPE, detail::layout_traits::Layout2d>
{
   using Layout = detail::layout_traits::Layout2d;

   SAMRAI_INLINE ArrayView<2, const TYPE>(const TYPE* data, const hier::Box& box) :
      RAJA::View<const TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 2>{ {box.lower()[0], box.lower()[1]} },
            std::array<RAJA::Index_type, 2>{ {box.upper()[0], box.upper()[1]} },
            RAJA::as_array<RAJA::PERM_JI>::get())){}
};

template<class TYPE>
struct ArrayView<3, const TYPE> : public RAJA::View<const TYPE, detail::layout_traits::Layout3d>
{
   using Layout = detail::layout_traits::Layout3d;

   SAMRAI_INLINE ArrayView<3, const TYPE>(const TYPE* data, const hier::Box& box) :
      RAJA::View<const TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 3>{ {box.lower()[0], box.lower()[1], box.lower()[2]} },
            std::array<RAJA::Index_type, 3>{ {box.upper()[0], box.upper()[1], box.upper()[2]} },
            RAJA::as_array<RAJA::PERM_KJI>::get())){};
};

}
}

#endif // HAVE_RAJA

#endif
