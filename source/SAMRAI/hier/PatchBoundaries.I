/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Describes boundaries for a patch
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *                                                                       *
 * Default constructor                                                   *
 *                                                                       *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
PatchBoundaries::PatchBoundaries():
   d_dim(tbox::Dimension::getInvalidDimension()),
   d_array_of_bboxes(0)
{
}

/*
 *************************************************************************
 *                                                                       *
 * Array access operator                                                 *
 *                                                                       *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
tbox::Array<hier::BoundaryBox>&
PatchBoundaries::operator [] (unsigned int i)
{
   TBOX_ASSERT(i < d_dim.getValue());
   return d_array_of_bboxes[i];
}

/*
 *************************************************************************
 *                                                                       *
 * Array access operator                                                 *
 *                                                                       *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Array<hier::BoundaryBox>&
PatchBoundaries::operator [] (unsigned int i) const
{
   TBOX_ASSERT(i < d_dim.getValue());
   return d_array_of_bboxes[i];
}

/*
 *************************************************************************
 *                                                                       *
 * Get copy of arrays                                                    *
 *                                                                       *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
tbox::Array<tbox::Array<hier::BoundaryBox> >
PatchBoundaries::getArrays()
{
   return d_array_of_bboxes;
}

/*
 *************************************************************************
 *                                                                       *
 * Get copy of arrays                                                    *
 *                                                                       *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Array<tbox::Array<hier::BoundaryBox> >
PatchBoundaries::getArrays() const
{
   return d_array_of_bboxes;
}

}
}
