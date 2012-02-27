/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Transformation class
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
Transformation::RotationIdentifier
Transformation::getRotation() const
{
   return d_rotation;
}

SAMRAI_INLINE_KEYWORD
const IntVector&
Transformation::getOffset() const
{
   return d_offset;
}

SAMRAI_INLINE_KEYWORD
Transformation&
Transformation::operator = (
   const Transformation& rhs)
{
   d_rotation = rhs.d_rotation;
   d_offset = rhs.d_offset;
   d_begin_block = rhs.d_begin_block;
   d_end_block = rhs.d_end_block;

   return *this;
}

/*
 *************************************************************************
 *
 * Rotate an index around the origin.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
Transformation::rotateIndex(
   Index& index,
   const RotationIdentifier rotation)
{
   rotateIndex(&index[0], index.getDim(), rotation);
}

}
}
