/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A list of boxes with basic domain calculus operations 
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxList::BoxList():
   tbox::List<Box>(),
   d_dim(tbox::Dimension::getInvalidDimension())
{
}

SAMRAI_INLINE_KEYWORD
BoxList::BoxList(const tbox::Dimension& dim):
   tbox::List<Box>(),
   d_dim(dim)
{
}

SAMRAI_INLINE_KEYWORD
int BoxList::getNumberOfBoxes() const
{
   return getNumberOfItems();
}

SAMRAI_INLINE_KEYWORD
int BoxList::size() const
{
   return getNumberOfItems();
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& BoxList::getDim() const
{
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
void BoxList::unionBoxes(
   const Box& box)
{
   if (!box.empty()) {
      appendItem(box);
   }
}

SAMRAI_INLINE_KEYWORD
void BoxList::unionBoxes(
   const BoxList& boxes)
{
   copyItems(boxes);
}

SAMRAI_INLINE_KEYWORD
BoxList::~BoxList()
{
}

}
}
