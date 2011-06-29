/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for node centered patch data types 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
NodeIterator::NodeIterator(
   const hier::Box& box):
   d_index(box.lower(), hier::IntVector::getZero(box.getDim())),
   d_box(NodeGeometry::toNodeBox(box))
{
}

SAMRAI_INLINE_KEYWORD
NodeIterator::NodeIterator(
   const NodeIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

SAMRAI_INLINE_KEYWORD
NodeIterator&
NodeIterator::operator = (
   const NodeIterator& iter)
{
   d_index = iter.d_index;
   d_box = iter.d_box;
   return *this;
}

SAMRAI_INLINE_KEYWORD
NodeIterator::~NodeIterator()
{
}

SAMRAI_INLINE_KEYWORD
const NodeIndex& NodeIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const NodeIndex& NodeIterator::operator () () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
NodeIterator::operator bool () const
{
   bool retval = true;
   for (int i = 0; i < d_box.getDim().getValue(); i++) {
      if (d_index(i) > d_box.upper(i)) {
         retval = false;
         break;
      }
   }

   return retval;
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
SAMRAI_INLINE_KEYWORD
NodeIterator::operator const void* () const
{
   return NodeIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool NodeIterator::operator ! () const
{
   return !NodeIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
void NodeIterator::operator ++ (
   int)
{
   d_index(0)++;
   for (int i = 0; i < d_box.getDim().getValue() - 1; i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i + 1)++;
      } else {
         break;
      }
   }
}

SAMRAI_INLINE_KEYWORD
bool NodeIterator::operator == (
   const NodeIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool NodeIterator::operator != (
   const NodeIterator& iter) const
{
   return d_index != iter.d_index;
}

}
}
