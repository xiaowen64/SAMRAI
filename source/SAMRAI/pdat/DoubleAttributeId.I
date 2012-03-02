/**********************************************************************
*
* This file is part of the SAMRAI distribution.  For full copyright
* information, see COPYRIGHT and COPYING.LESSER.
*
* Copyright:     (c) 1997 - 2011 Lawrence Livermore National Security, LLC
* Description:   pdat
**********************************************************************/

namespace SAMRAI {
namespace pdat {

/**********************************************************************
 * DoubleAttributeId equality
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
bool
DoubleAttributeId::operator == (
   const DoubleAttributeId& rhs) const
{
   return d_val == rhs.d_val;
}

/**********************************************************************
 * DoubleAttributeId inequality
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
bool
DoubleAttributeId::operator != (
   const DoubleAttributeId& rhs) const
{
   return !this->operator == (rhs);
}

/**********************************************************************
*
**********************************************************************/
SAMRAI_INLINE_KEYWORD
int
DoubleAttributeId::operator () () const
{
   return d_val;
}

/**********************************************************************
* Assignment
**********************************************************************/
SAMRAI_INLINE_KEYWORD
DoubleAttributeId&
DoubleAttributeId::operator = (
   const DoubleAttributeId& rhs)
{
   if (this != &rhs) {
      d_val = rhs.d_val;
   }
   return *this;
}

}
}
