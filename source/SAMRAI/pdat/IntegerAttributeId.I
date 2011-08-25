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
*
**********************************************************************/
SAMRAI_INLINE_KEYWORD
int
IntegerAttributeId::operator () () const
{
   return d_val;
}

/**********************************************************************
 * IntegerAttributesId ctor
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
IntegerAttributeId::IntegerAttributeId(
   int value):
   d_val(value)
{
}

/**********************************************************************
 * IntegerAttributesId copy ctor
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
IntegerAttributeId::IntegerAttributeId(
   const IntegerAttributeId& other):
   d_val(other.d_val)
{
}

/**********************************************************************
 * d'tor
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
IntegerAttributeId::~IntegerAttributeId()
{
}

/**********************************************************************
 * IntegerAttributesId equality
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
bool
IntegerAttributeId::operator == (
   const IntegerAttributeId& rhs) const
{
   return d_val == rhs.d_val;
}

/**********************************************************************
 * IntegerAttributesId inequality
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
bool
IntegerAttributeId::operator != (
   const IntegerAttributeId& rhs) const
{
   return !this->operator == (rhs);
}

/**********************************************************************
 * IntegerAttributesId c'tor (private)
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
IntegerAttributeId::IntegerAttributeId():
   d_val(-1)
{
}

/**********************************************************************
 * assignment
 *********************************************************************/
SAMRAI_INLINE_KEYWORD
IntegerAttributeId&
IntegerAttributeId::operator = (
   const IntegerAttributeId& rhs)
{
   if (this != &rhs) {
      d_val = rhs.d_val;
   }
   return *this;
}

}
}
