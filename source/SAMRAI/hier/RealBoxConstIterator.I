/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator over real Boxes in a BoxContainer.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
bool
RealBoxConstIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_ni != d_mapped_boxes->end() &&
          !d_ni->isPeriodicImage();
}

SAMRAI_INLINE_KEYWORD
RealBoxConstIterator&
RealBoxConstIterator::operator = (
   const RealBoxConstIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_ni = r.d_ni;
   return *this;
}

SAMRAI_INLINE_KEYWORD
const Box&
RealBoxConstIterator::operator * () const
{
   return *d_ni;
}

SAMRAI_INLINE_KEYWORD
const Box*
RealBoxConstIterator::operator -> () const
{
   return &(*d_ni);
}

SAMRAI_INLINE_KEYWORD
bool
RealBoxConstIterator::operator == (
   const RealBoxConstIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes && d_ni == r.d_ni;
}

SAMRAI_INLINE_KEYWORD
bool
RealBoxConstIterator::operator != (
   const RealBoxConstIterator& r) const
{
   return d_mapped_boxes != r.d_mapped_boxes || d_ni != r.d_ni;
}

}
}
